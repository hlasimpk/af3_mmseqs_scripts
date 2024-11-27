#!/usr/bin/env python

import os
import json

from af3_script_utils import (
    custom_template_argpase_util,
    get_custom_template,
)


def run_custom_template(
    input_json,
    target_id,
    custom_template,
    custom_template_chain,
    output_json=None,
    to_file=True,
):
    af3_json = json.load(open(input_json))
    if not os.path.exists(custom_template):
        msg = f"Custom template file {custom_template} not found"
        raise FileNotFoundError(msg)

    for sequence in af3_json["sequences"]:
        if "protein" not in sequence:
            continue

        sequence = get_custom_template(
            sequence,
            target_id,
            custom_template,
            custom_template_chain,
        )

    if to_file:
        if not output_json:
            output_json = input_json

        with open(output_json, "w") as f:
            json.dump(af3_json, f)

    return af3_json


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Add custom template to alphafold input JSON"
    )

    parser.add_argument("--input_json", help="Input alphafold3 json file")
    parser.add_argument("--output_json", help="Output alphafold3 json file")
    parser = custom_template_argpase_util(parser)

    args = parser.parse_args()

    run_custom_template(
        args.input_json,
        args.target_id,
        args.custom_template,
        args.custom_template_chain,
        output_json=args.output_json,
        to_file=True,
    )
