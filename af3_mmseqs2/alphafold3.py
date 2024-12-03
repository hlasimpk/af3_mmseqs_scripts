# start by finding the directory where the alphafold3.py script is located

from af3_mmseqs2.add_custom_template import custom_template_argpase_util
from af3_mmseqs2.add_mmseqs_msa import mmseqs2_argparse_util, add_msa_to_json
from af3_mmseqs2.af3_script_utils import setup_logger
import configparser
import json
from pathlib import Path
import subprocess
import sys

logger = setup_logger()

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
    
    logger.info("Running Alphafold3")
    with subprocess.Popen(
        cmd, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE
    ) as p:
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            logger.error(stderr.decode())
            raise subprocess.CalledProcessError(p.returncode, cmd, stderr)
    logger.info("Alphafold3 run complete")
    logger.info("Output files are in %s", output_dir)


def af3_argparse_main(parser):
    parser.add_argument("input_json", help="Input sequence file")

    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--output_json", help="Output json file")
    # make the vartible saved as database_dir
    parser.add_argument(
        "--database",
        help="The Database directory for the generation of the MSA.",
        dest="database_dir",
        default=None,
    )
    parser.add_argument(
        "--mmseqs2",
        help="Use MMseqs2 for MSA",
        action="store_true",
    )

    parser.add_argument(
        "--model_params",
        help="The directory containing the model parameters",
        default=None,
    )

    mmseqs2_argparse_util(parser)
    custom_template_argpase_util(parser)

    return parser


def main():
    """Run AlphaFold3"""
    import argparse
    parser = argparse.ArgumentParser(description="Run AlphaFold3")

    # Load defaults from config file
    defaults = {}
    config_file = Path(__file__).parent.joinpath("..", "config.ini")
    config = configparser.SafeConfigParser()

    if config_file.exists():
        config.read(str(config_file))
        defaults.update(dict(config.items("Databases")))

    parser = af3_argparse_main(parser)
    parser.set_defaults(**defaults)
    args = parser.parse_args()

    updated_config = False
    if args.model_params != defaults["model_params"]:
        config.set("Databases", "model_params", args.model_params)
        updated_config = True
    if args.database_dir != defaults["database_dir"]:
        config.set("Databases", "database_dir", args.database_dir)
        updated_config = True
    if updated_config:
        with open(config_file, "w") as f:
            config.write(f)

    if not args.database_dir or not Path(args.database_dir).exists():
        logger.error(f"Database directory not found: {args.database_dir}")
        sys.exit(1)
    elif not args.model_params or not Path(args.model_params).exists():
        logger.error(f"Model parameters directory not found: {args.model_params}")
        sys.exit(1)

    with open(args.input_json, "r") as f:
        af3_json = json.load(f)

    if args.mmseqs2:
        af3_json = add_msa_to_json(
            input_json=args.input_json,
            templates=args.templates,
            num_templates=args.num_templates,
            custom_template=args.custom_template,
            custom_template_chain=args.custom_template_chain,
            target_id=args.target_id,
            af3_json=af3_json,
            output_json=args.output_json,
            to_file=True,
        )

        output_json = (
            args.input_json.replace(".json", "_mmseqs.json")
            if args.output_json is None
            else args.output_json
        )
    else:
        output_json = args.input_json
    run_alphafold3(
        input_json=output_json,
        output_dir=args.output_dir,
        model_params=args.model_params,
        database_dir=args.database_dir,
    )


if __name__ == "__main__":
    main()
