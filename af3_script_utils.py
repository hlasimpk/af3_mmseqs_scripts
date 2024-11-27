from Bio.PDB import MMCIFParser, MMCIFIO
from Bio import pairwise2
from io import StringIO

import os
import time


def check_chains(mmcif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("template", mmcif_file)
    chains = []
    for model in structure:
        for chain in model:
            chains.append(chain.id)
    return chains


def extract_sequence_from_mmcif(mmcif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("template", mmcif_file)
    sequence = ""
    for model in structure:
        for chain in model:  # Assuming one chain only
            for residue in chain:
                if residue.id[0] == " ":  # Exclude heteroatoms
                    sequence += residue.resname[
                        0
                    ]  # Simplified to take the first letter
    return sequence


def align_and_map(query_seq, template_seq):
    # Perform pairwise alignment
    alignments = pairwise2.align.globalxx(query_seq, template_seq)
    alignment = alignments[0]  # Take the best alignment
    query_aligned, template_aligned, _, _, _ = alignment

    # Map indices
    query_indices = []
    template_indices = []

    query_pos, template_pos = -1, -1
    for q_char, t_char in zip(query_aligned, template_aligned):
        if q_char != "-":  # Not a gap in query
            query_pos += 1
            if t_char != "-":  # Not a gap in template
                query_indices.append(query_pos)
                # 1-based index for template
                template_indices.append(template_pos + 1)
        elif t_char != "-":  # Increment template position for gaps in query
            template_pos += 1

    return query_indices, template_indices


def get_mmcif(
    cif,
    pdb_id,
    chain_id,
    start,
    end,
    tmpdir=None,
):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, cif)

    # Extract release date from the CIF file
    mmcif_dict = parser._mmcif_dict
    headers_to_keep = [
        "_entry.id",
        "_entry.title",
        "_entry.deposition_date",
        "_pdbx_audit_revision_history.revision_date",
    ]
    filtered_metadata = {
        key: mmcif_dict[key] for key in headers_to_keep if key in mmcif_dict
    }

    # Make metadata if missing
    if "_pdbx_audit_revision_history.revision_date" not in filtered_metadata:
        filtered_metadata["_pdbx_audit_revision_history.revision_date"] = time.strftime(
            "%Y-%m-%d"
        )

    for model in structure:
        chain_to_del = []
        for chain in model:
            if chain.id != chain_id:
                chain_to_del.append(chain.id)
                continue

        for unwanted_chain in chain_to_del:
            model.detach_child(unwanted_chain)

        for chain in model:
            res_to_del = []
            for i, res in enumerate(chain):
                rel_pos = i + 1
                if rel_pos < start or rel_pos > end or res.id[0] != " ":
                    res_to_del.append(res)

            for res in res_to_del:
                chain.detach_child(res.id)

    # Save the filtered structure to a new CIF file
    io = MMCIFIO()
    io.set_structure(structure)
    filtered_output = (
        f"{pdb_id}_{chain_id}.cif"
        if tmpdir is None
        else f"{tmpdir}/{pdb_id}_{chain_id}.cif"
    )
    io.save(filtered_output)

    # Parse the filtered structure to get the modified MMCIF with no metadata
    structure = parser.get_structure(pdb_id, filtered_output)
    mmcif_dict = parser._mmcif_dict

    # Add the filtered metadata to the MMCIF dictionary
    mmcif_dict.update(filtered_metadata)

    # Save the modified MMCIF with wanted metadata to a string
    string_io = StringIO()
    io.set_dict(mmcif_dict)
    io.save(string_io)

    os.unlink(filtered_output)

    return string_io.getvalue()


def mmseqs2_argparse_util(parser):
    parser.add_argument(
        "--templates", action="store_true", help="Include templates in the output json"
    )
    parser.add_argument(
        "--num_templates",
        type=int,
        default=20,
        help="Number of templates to include in the output json",
    )
    return parser


def custom_template_argpase_util(parser):
    parser.add_argument("--target_id", help="Target id relating to the custom template")
    parser.add_argument(
        "--custom_template", help="Custom template to include in the output json"
    )
    parser.add_argument(
        "--custom_template_chain",
        help="Custom template chain to include in the output json",
    )

    return parser


# runs for each sequence in the input json
def get_custom_template(
    sequence,
    target_id,
    custom_template,
    custom_template_chain,
):
    # add code here to run the custom template but move the res to output
    # Keep existing templates if they're present
    if "templates" not in sequence["protein"]:
        templates = []
    else:
        templates = sequence["protein"]["templates"]

    input_sequence = sequence["protein"]["sequence"]
    seq_id = sequence["protein"]["id"]
    if isinstance(seq_id, list):
        if target_id and target_id not in seq_id:
            return sequence
    if isinstance(seq_id, str):
        if target_id and target_id != seq_id:
            return sequence

    if not os.path.exists(custom_template):
        msg = f"Custom template file {custom_template} not found"
        raise FileNotFoundError(msg)

    chain_info = check_chains(custom_template)
    if len(chain_info) != 1 and not custom_template_chain:
        msg = f"Custom template file {custom_template} contains \
{len(chain_info)} chains. Please specify the chain to use with --custom_template_chain"
        raise ValueError(msg)

    if custom_template_chain and custom_template_chain not in chain_info:
        msg = f"Custom template file {custom_template} does not \
contain chain {custom_template_chain}"
        raise ValueError(msg)

    if not custom_template_chain:
        custom_template_chain = chain_info[0]

    template = {}
    cif_str = get_mmcif(
        custom_template,
        "custom",
        custom_template_chain,
        1,
        len(input_sequence),
    )

    template["mmcif"] = cif_str
    template_seq = extract_sequence_from_mmcif(StringIO(cif_str))
    query_indices, template_indices = align_and_map(input_sequence, template_seq)

    template["queryIndices"] = query_indices
    template["templateIndices"] = template_indices

    # Add the custom template to the start of the templates list
    templates.insert(0, template)

    # Add template to the json
    sequence["protein"]["templates"] = templates

    # Save the output json
    return sequence
