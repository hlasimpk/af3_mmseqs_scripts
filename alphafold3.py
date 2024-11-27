from add_custom_template import custom_template_argpase_util
from add_mmseqs_msa import mmseqs2_argparse_util, add_msa_to_json
import json
from pathlib import Path
import subprocess


def run_alphafold3(
    input_json: str | Path,
    output_dir: str | Path,
    model_params: str | Path,
    database_dir: str | Path,
) -> None:
    input_json = Path(input_json)
    output_dir = Path(output_dir)
    cmd = rf"""
    docker run -it \
    --volume {input_json.parent.resolve()}:/root/af_input \
    --volume {output_dir.resolve()}:/root/af_output \
    --volume {model_params}:/root/models \
    --volume {database_dir}:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input/{input_json.name} \
    --model_dir=/root/models \
    --output_dir=/root/af_output
    """

    with subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as p:
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            print(stderr.decode())
            print(stdout.decode())
            raise subprocess.CalledProcessError(p.returncode, cmd, stderr)

    print(stdout.decode())
    print(stderr.decode())

    print("Alphafold3 run complete")
    print("Output files are in", output_dir)


def af3_argparse_main(parser):
    parser.add_argument("input_json", help="Input sequence file", nargs="+")

    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--output_json", help="Output json file", nargs="+")
    # make the vartible saved as database_dir
    parser.add_argument(
        "--database",
        help="The Database directory for the generation of the MSA.",
        required=True,
        dest="database_dir",
    )
    parser.add_argument(
        "--mmseqs2",
        help="Use MMseqs2 for MSA",
        action="store_true",
    )

    parser.add_argument(
        "--model_params",
        help="The directory containing the model parameters",
        required=True,
    )

    mmseqs2_argparse_util(parser)
    custom_template_argpase_util(parser)

    return parser


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run AlphaFold3")

    parser = af3_argparse_main(parser)

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
        with open(json_file, "r") as f:
            af3_json = json.load(f)

            if args.mmseqs2:
                af3_json = add_msa_to_json(
                    input_json=json_file,
                    templates=args.templates,
                    num_templates=args.num_templates,
                    custom_template=args.custom_template,
                    custom_template_chain=args.custom_template_chain,
                    target_id=args.target_id,
                    af3_json=af3_json,
                    output_json=args.output_json[i],
                    to_file=True,
                )

                run_json = (
                    json_file.replace(".json", "_mmseqs.json")
                    if args.output_json[i] is None
                    else args.output_json[i]
                )
            else:
                run_json = json_file

        run_alphafold3(
            input_json=run_json,
            output_dir=args.output_dir,
            model_params=args.model_params,
            database_dir=args.database_dir,
        )
