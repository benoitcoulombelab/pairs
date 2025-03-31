import argparse
import glob
import json
import os
import sys
from typing import TextIO

import tqdm

RANKING_METRICS = ["interface", "pitms", "plddts", "ptms"]
RANKING_METRICS_JSON = {
  "interface": "interface score",
  "pitms": "pitms",
  "plddts": "plddts",
  "ptms": "ptms"
}
RANKING_FILES = ["ranking_debug.json", "ranking_all_*.json",
                 "ranking_model_*.json"]


def dir_path(string: str):
  if not string or os.path.isdir(string):
    return string
  else:
    raise NotADirectoryError(string)


def main(argv: list[str] = None):
  parser = argparse.ArgumentParser(
      description="List files from AlphaFold's output that need to be archived.")
  parser.add_argument('input', nargs='?', type=dir_path, default="",
                      help="Directory containing one or many AlphaFold's output directories  "
                           "(default: current directory)")
  parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout,
                      help="Output file  (default: standard output)")
  parser.add_argument('-p', '--progress', action="store_true", default=False,
                      help="Show progress bar")
  parser.add_argument('-m', '--metric', choices=RANKING_METRICS,
                      default=RANKING_METRICS[0],
                      help="Metric to use to choose best model - only useful for AF2Complex runs "
                           "as the metric is always 'plddts' for AlphaFold runs  (default: %(default)s)")
  parser.add_argument('--all-pdb', action="store_true", default=False,
                      help="Archive all PDB files")
  parser.add_argument('--pkl', action=argparse.BooleanOptionalAction,
                      default=True,
                      help="Archive PKL files of the best model")
  parser.add_argument('--all-pkl', action="store_true", default=False,
                      help="Archive all models' PKL files")

  args = parser.parse_args(argv)

  list_files(input_dir=args.input, output_file=args.output,
             progress=args.progress, metric=args.metric,
             all_pdb=args.all_pdb, best_pkl=args.pkl, all_pkl=args.all_pkl)


def list_files(input_dir: str, output_file: TextIO, progress: bool = False,
    metric: str = RANKING_METRICS[0],
    all_pdb: bool = False,
    best_pkl: bool = True, all_pkl: bool = False):
  """
  List files from AlphaFold's output that need to be archived.

  :param input_dir: directory containing one or many AlphaFold's output directories
  :param output_file: output file
  :param progress: show progress bar
  :param metric: metric used to choose best model - only used for AF2Complex rankings
  :param all_pdb: archive all PDB files
  :param best_pkl: archive PKL files of the best model
  :param all_pkl: archive all models' PKL files
  """
  ranking_files = []
  [ranking_files.extend(
      glob.glob(os.path.join(input_dir, f"**/{ranking_file}"))) for ranking_file
    in RANKING_FILES]
  directories = [os.path.dirname(ranking_file) for ranking_file in
                 ranking_files]
  directories = list(set(directories))
  directories.sort()
  for directory in (tqdm.tqdm(directories) if progress else directories):
    ranking_files = []
    [ranking_files.extend(glob.glob(os.path.join(directory, ranking_file))) for
     ranking_file in RANKING_FILES]
    ranking_files.sort()
    best_model = get_best_model(ranking_files, metric)
    files_to_archive = glob.glob(os.path.join(directory, "*.json"))
    if all_pdb:
      files_to_archive.extend(glob.glob(os.path.join(directory, "*.pdb")))
    else:
      files_to_archive.extend(
          glob.glob(os.path.join(directory, f"*{best_model}*.pdb")))
      files_to_archive.extend(
          glob.glob(os.path.join(directory, f"ranked_0.pdb")))
    if all_pkl:
      files_to_archive.extend(
          glob.glob(os.path.join(directory, "*model_*.pkl")))
    elif best_pkl:
      files_to_archive.extend(
          glob.glob(os.path.join(directory, f"*{best_model}*.pkl")))
    files_to_archive = list(set(files_to_archive))
    files_to_archive.sort()
    for file in files_to_archive:
      output_file.write(file)
      output_file.write("\n")


def get_best_model(ranking_files: [str],
    metric: str = RANKING_METRICS[0]) -> str:
  """
  Returns best model found by AlphaFold or AF2Complex in ranking file.

  :param ranking_files: ranking files in JSON format - usually named 'ranking_debug.json',
                        'ranking_all_*.json' or 'ranking_model_*.json'
  :param metric: metric used to choose best model - only used for AF2Complex rankings
  :return: best model found by AlphaFold or AF2Complex in ranking file
  """
  if metric not in RANKING_METRICS_JSON:
    raise AssertionError(
        f"metric {metric} not found in RANKING_METRICS ({RANKING_METRICS.keys()})")
  ranking_json = RANKING_METRICS_JSON[metric]
  rankings = []
  for ranking_file in ranking_files:
    with open(ranking_file, 'r') as ranking_in:
      rankings.append(json.load(ranking_in))
  if len(rankings) == 1 and ranking_json not in rankings[0]:
    # Assume AlphaFold ranking
    return rankings[0]["order"][0]
  else:
    # Assume AF2Complex ranking
    best_score = -1.0
    best_model = None
    for ranking in rankings:
      all_keys = list(ranking[ranking_json].keys())
      keys = [key for key in all_keys if "_recycled_" not in key]
      for key in keys:
        score = ranking[ranking_json][key]
        if score > best_score:
          best_score = score
          best_model = key
    return best_model


if __name__ == '__main__':
  main()
