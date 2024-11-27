#!/usr/bin/env python

import os
import time
import tarfile
import random
import logging
import json
from io import StringIO
from typing import Tuple, List

import requests
from tqdm.autonotebook import tqdm

import tempfile
from af3_script_utils import (
    custom_template_argpase_util,
    mmseqs2_argparse_util,
    align_and_map,
    extract_sequence_from_mmcif,
    get_mmcif,
    get_custom_template,
)


logger = logging.getLogger(__name__)

TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"
)


class MMseqs2Exception(Exception):
    def __init__(self):
        msg = "MMseqs2 API is giving errors. Please confirm your input is a valid \
protein sequence. If error persists, please try again an hour later."
        super().__init__(msg)


def add_msa_to_json(
    input_json,
    templates,
    num_templates,
    custom_template,
    custom_template_chain,
    target_id,
    af3_json=None,
    output_json=None,
    to_file=True,
):
    if af3_json is None:
        with open(input_json, "r") as f:
            af3_json = json.load(f)

    with tempfile.TemporaryDirectory() as tmpdir:
        for sequence in af3_json["sequences"]:
            if "protein" in sequence:
                input_sequence = sequence["protein"]["sequence"]

                # Run MMseqs2 to get unpaired MSA
                if templates:
                    a3m_lines, templates = run_mmseqs(
                        input_sequence,
                        tmpdir,
                        use_templates=True,
                        num_templates=num_templates,
                    )

                else:
                    a3m_lines = run_mmseqs(input_sequence, tmpdir, use_templates=False)
                    templates = []

                if custom_template:
                    if not os.path.exists(custom_template):
                        msg = f"Custom template file {custom_template} not found"
                        raise FileNotFoundError(msg)
                    # Can only add templates to protein sequences, so check if there
                    # are multiple protein sequences in the input json
                    if (
                        len([x for x in af3_json["sequences"] if "protein" in x.keys()])
                        > 1
                        and not target_id
                    ):
                        msg = "Multiple sequences found in input json. Please specify \
target id so that custom template can be added to the correct sequence"
                        raise ValueError(msg)

                    sequence = get_custom_template(
                        sequence,
                        target_id,
                        custom_template,
                        custom_template_chain,
                    )

                # Add unpaired MSA to the json
                sequence["protein"]["unpairedMsa"] = a3m_lines[0]
                sequence["protein"]["pairedMsa"] = ""
                sequence["protein"]["templates"] = templates
    if to_file:
        if output_json:
            with open(output_json, "w") as f:
                json.dump(af3_json, f)
        else:
            output_json = input_json.replace(".json", "_mmseqs.json")
            with open(output_json, "w") as f:
                json.dump(af3_json, f)

    return af3_json


def run_mmseqs(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    host_url="https://a3m.mmseqs.com",
    num_templates=20,
) -> Tuple[List[str], List[str]]:
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        res = requests.post(
            f"{host_url}/{submission_endpoint}", data={"q": query, "mode": mode}
        )
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        res = requests.get(f"{host_url}/ticket/{ID}")
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        res = requests.get(f"{host_url}/result/download/{ID}")
        with open(path, "wb") as out:
            out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        mode = ""
        use_templates = False
        use_env = False

    # define path
    path = prefix
    if not os.path.isdir(path):
        os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = list(set(seqs))
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise MMseqs2Exception()

                if out["status"] == "MAINTENANCE":
                    raise MMseqs2Exception()

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise MMseqs2Exception()

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

    tested_pdbs = []
    templates = []
    if use_templates:
        count = 0
        for line in open(f"{path}/pdb70.m8", "r"):
            template = {}
            if count < num_templates:
                p = line.rstrip().split()
                M, pdb, qid, alilen, tstart, tend = (
                    p[0],
                    p[1],
                    p[2],
                    p[3],
                    p[8],
                    p[9],
                )
                coverage = float(alilen) / len(x)
                pdb_id = pdb.split("_")[0]

                # Use the same template filters as AF3 and only use 1 template per PDB
                if (
                    float(qid) == 1.0
                    and coverage >= 0.95
                    or coverage < 0.1
                    or pdb_id in tested_pdbs
                ):
                    continue

                pdb_id = pdb.split("_")[0]
                cif_str = fetch_mmcif(
                    pdb_id, pdb.split("_")[1], int(tstart), int(tend), prefix
                )
                template["mmcif"] = cif_str

                template_seq = extract_sequence_from_mmcif(StringIO(cif_str))
                query_indices, template_indices = align_and_map(x, template_seq)

                template["queryIndices"] = query_indices
                template["templateIndices"] = template_indices
                templates.append(template)
                tested_pdbs.append(pdb_id)
                count += 1

    return (a3m_lines, templates) if use_templates else a3m_lines


def fetch_mmcif(
    pdb_id,
    chain_id,
    start,
    end,
    tmpdir,
):
    pdb_id = pdb_id.lower()
    url_base = "http://www.ebi.ac.uk/pdbe-srv/view/files/"
    url = url_base + pdb_id + ".cif"
    response = requests.get(url)
    text = response.text

    output = tmpdir + pdb_id + ".cif"
    with open(output, "w") as f:
        f.write(text)

    return get_mmcif(output, pdb_id, chain_id, start, end, tmpdir)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Add MMseqs2 unpaired MSA to AlphaFold3 json"
    )
    parser.add_argument("--input_json", help="Input alphafold3 json file", nargs="+")
    parser.add_argument("--output_json", help="Output alphafold3 json file", nargs="+")

    parser = mmseqs2_argparse_util(parser)
    parser = custom_template_argpase_util(parser)

    args = parser.parse_args()
    if args.output_json and len(args.input_json) != len(args.output_json):
        msg = "If output_json is specified, the number of output json files must \
match the number of input json files"
        raise ValueError(msg)

    if len(args.input_json) > 1 and args.custom_template:
        msg = "Multiple input json files found. This is not supported with custom \
template. Please run custom template separately for each input json file"
        raise ValueError(msg)

    if not args.output_json:
        args.output_json = [None] * len(args.input_json)

    for i, json_file in enumerate(args.input_json):
        add_msa_to_json(
            json_file,
            args.templates,
            args.num_templates,
            args.custom_template,
            args.custom_template_chain,
            args.target_id,
            output_json=args.output_json[i],
            to_file=True,
        )
