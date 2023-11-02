import argparse
import os
import re
import sys
import glob
import json
from typing import TextIO

import tqdm


RELAXATION_FILENAMES = "relaxed_model_*.pdb"
RELAXATION_FILENAME_TO_MODEL_REGEX = r"relaxed_(model_[\w]+)\.pdb"
RELAXATION_JSON = "relax_metrics.json"
ALPHAFOLD_FILES_TO_ARCHIVE = ["ranked_0.pdb", "ranking_debug.json", "relax_metrics.json", "timings.json"]
ALPHAFOLD_FILES_TO_ARCHIVE.extend([f"unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb"
                                   for i in range(1, 6) for j in range(0, 5)])


def dir_path(string: str):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="List files from AlphaFold's output that need to be archived.")
    parser.add_argument('input', nargs='?', type=dir_path, default="",
                        help="Directory containing one or many AlphaFold's output directories  "
                             "(default: current directory)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file  (default: standard output)")
    parser.add_argument('-p', '--progress', action="store_true", default=False,
                        help="Show progress bar")

    args = parser.parse_args(argv)

    list_files(input_dir=args.input, output_file=args.output, progress=args.progress)


def list_files(input_dir: str, output_file: TextIO, progress: bool = False):
    """
    List files from AlphaFold's output that need to be archived.

    :param input_dir: directory containing one or many AlphaFold's output directories
    :param output_file: output file
    :param progress: show progress bar
    """
    ranked_0_files = glob.glob(os.path.join(input_dir, "**/ranked_0.pdb"))
    for ranked_0_file in (tqdm.tqdm(ranked_0_files) if progress else ranked_0_files):
        files = alphafold_archive_files(os.path.dirname(ranked_0_file))
        for file in files:
            output_file.write(file)
            output_file.write("\n")


def alphafold_archive_files(input_dir: str) -> list[str]:
    """
    List files from AlphaFold's output directory that need to be archived.

    :param input_dir: AlphaFold's output directory
    :return: all files to be archived
    """
    files_to_archive = []
    for file in ALPHAFOLD_FILES_TO_ARCHIVE:
        file = os.path.join(input_dir, file)
        if os.path.isfile(file):
            files_to_archive.append(file)
        else:
            print(f"File {os.path.basename(file)} is missing from directory {input_dir}", file=sys.stderr)
    relaxation_models_to_archive = relaxation_models(input_dir)
    for relaxation_model in relaxation_models_to_archive:
        for file in [f"relaxed_{relaxation_model}.pdb", f"result_{relaxation_model}.pkl"]:
            file = os.path.join(input_dir, file)
            if os.path.isfile(file):
                files_to_archive.append(file)
            else:
                print(f"File {os.path.basename(file)} is missing from directory {input_dir}", file=sys.stderr)
    return files_to_archive


def relaxation_models(input_dir: str) -> list[str]:
    """
    Returns all models that were relaxed.
    When possible, models are obtained using relax_metrics.json file.
    Otherwise, models are obtained based on filename.

    :param input_dir: AlphaFold's output directory
    :return: all models that were relaxed
    """
    relaxation_json = os.path.join(input_dir, RELAXATION_JSON)
    if os.path.isfile(relaxation_json):
        with open(relaxation_json, 'r') as relaxation_in:
            relaxation = json.load(relaxation_in)
        return relaxation.keys()
    relaxation_files = glob.glob(os.path.join(input_dir, RELAXATION_FILENAMES))
    return [re.match(RELAXATION_FILENAME_TO_MODEL_REGEX, os.path.basename(relaxation_file))[1]
            for relaxation_file in relaxation_files]


if __name__ == '__main__':
    main()
