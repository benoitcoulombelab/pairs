import os
import shutil
import sys
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import pytest

from pairs import ListFiles


@pytest.fixture
def mock_testclass():
  _list_files = ListFiles.list_files
  _get_best_model = ListFiles.get_best_model
  yield
  ListFiles.list_files = _list_files
  ListFiles.get_best_model = _get_best_model


def create_alphafold3_files(alphafold_output, name):
  create_files = [f"{name}_confidences.json", f"{name}_data.json",
                  f"{name}_model.cif",
                  f"{name}_summary_confidences.json", "ranking_scores.csv",
                  "TERMS_OF_USE.md"]
  [open(os.path.join(alphafold_output, file), 'w') for file in create_files]
  sample_folders = ["seed-1_sample-0", "seed-1_sample-1", "seed-1_sample-2",
                    "seed-1_sample-3", "seed-1_sample-4"]
  create_sample_files = ["confidences.json", "model.cif",
                         "summary_confidences.json"]
  for sample_folder in sample_folders:
    os.mkdir(os.path.join(alphafold_output, sample_folder))
    [open(os.path.join(alphafold_output, sample_folder, file), 'w') for file in
     create_sample_files]


def create_alphafold_files(alphafold_output):
  os.mkdir(os.path.join(alphafold_output, "msas"))
  os.mkdir(os.path.join(alphafold_output, "msas/A"))
  os.mkdir(os.path.join(alphafold_output, "msas/B"))
  create_files = ["features.pkl", "ranking_debug.json", "relax_metrics.json",
                  "timings.json",
                  "msas/chain_id_map.json"]
  create_files.extend([f"ranked_{i}.pdb" for i in range(0, 25)])
  create_files.extend(
      [f"unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" for i in range(1, 6) for
       j in range(0, 5)])
  create_files.extend(
      [f"relaxed_model_{i}_multimer_v3_pred_{j}.pdb" for i in range(1, 6) for j
       in range(0, 5)])
  create_files.extend(
      [f"result_model_{i}_multimer_v3_pred_{j}.pkl" for i in range(1, 6) for j
       in range(0, 5)])
  create_msas_files = ["mgnify_hits.sto", "small_bfd_hits.sto",
                       "uniref90_hits.sto",
                       "pdb_hits.sto", "uniprot_hits.sto"]
  [open(os.path.join(alphafold_output, file), 'w') for file in create_files]
  [open(os.path.join(alphafold_output, "msas/A", file), 'w') for file in
   create_msas_files]
  [open(os.path.join(alphafold_output, "msas/B", file), 'w') for file in
   create_msas_files]


def create_af2complex_files(af2complex_output, run_ids):
  os.mkdir(os.path.join(af2complex_output, "checkpoint"))
  os.mkdir(os.path.join(af2complex_output, "recycled"))
  if len(run_ids) == 1:
    create_files = [f"ranking_all_{run_ids[0]}.json", "relax_metrics.json"]
    create_files.extend(
        [f"model_{i}_multimer_v3_p1_{run_ids[0]}.pdb" for i in range(1, 6)])
    create_files.extend(
        [f"relaxed_model_{i}_multimer_v3_p1_{run_ids[0]}.pdb" for i in
         range(1, 6)])
    create_files.extend(
        [f"model_{i}_multimer_v3_p1_{run_ids[0]}.pkl" for i in range(1, 6)])
    create_files.extend(
        [f"recycled/model_{i}_multimer_v3_p1_{run_ids[0]}_recycled_0{j}.pdb"
         for i in range(1, 6) for j in range(0, 5)])
  else:
    create_files = ["relax_metrics.json"]
    create_files.extend(
        [f"ranking_model_{i}_multimer_v3_p1_{run_ids[i - 1]}.json" for i in
         range(1, 6)])
    create_files.extend(
        [f"model_{i}_multimer_v3_p1_{run_ids[i - 1]}.pdb" for i in range(1, 6)])
    create_files.extend(
        [f"relaxed_model_{i}_multimer_v3_p1_{run_ids[i - 1]}.pdb" for i in
         range(1, 6)])
    create_files.extend(
        [f"model_{i}_multimer_v3_p1_{run_ids[i - 1]}.pkl" for i in range(1, 6)])
    create_files.extend(
        [f"recycled/model_{i}_multimer_v3_p1_{run_ids[i - 1]}_recycled_0{j}.pdb"
         for i in range(1, 6) for j in range(0, 5)])
  create_files.extend(
      [f"checkpoint/model_{i}_multimer_v3_p1_False.pkl" for i in range(1, 6)])
  create_files.extend(
      [f"checkpoint/model_{i}_multimer_v3_p1_checkpoint.pkl" for i in
       range(1, 6)])
  [open(os.path.join(af2complex_output, file), 'w') for file in create_files]


def test_main(testdir, mock_testclass):
  ListFiles.list_files = MagicMock()
  ListFiles.main([])
  ListFiles.list_files.assert_called_once_with(
      input_dir="", output_file=ANY, progress=False, metric="interface",
      all_pdb=False, best_pkl=True, all_pkl=False)
  output_file = ListFiles.list_files.call_args.kwargs["output_file"]
  assert isinstance(output_file, TextIOWrapper)
  assert output_file.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
  input_dir = "alphafold"
  output_file = "output.txt"
  testdir.mkdir(input_dir)
  metric = "pitms"
  ListFiles.list_files = MagicMock()
  ListFiles.main(["-o", output_file, "-p", "-m", metric, input_dir])
  ListFiles.list_files.assert_called_once_with(
      input_dir=input_dir, output_file=ANY, progress=True, metric=metric,
      all_pdb=False, best_pkl=True, all_pkl=False)
  output_file_arg = ListFiles.list_files.call_args.kwargs["output_file"]
  assert isinstance(output_file_arg, TextIOWrapper)
  assert output_file_arg.mode in ["r+", "w"]
  assert output_file_arg.name == output_file


def test_main_long_parameters(testdir, mock_testclass):
  input_dir = "alphafold"
  output_file = "output.txt"
  testdir.mkdir(input_dir)
  metric = "pitms"
  ListFiles.list_files = MagicMock()
  ListFiles.main(
      ["--output", output_file, "--progress", "--metric", metric, "--all-pdb",
       "--no-pkl", "--all-pkl", input_dir])
  ListFiles.list_files.assert_called_once_with(
      input_dir=input_dir, output_file=ANY, progress=True, metric=metric,
      all_pdb=True, best_pkl=False, all_pkl=True)
  output_file_arg = ListFiles.list_files.call_args.kwargs["output_file"]
  assert isinstance(output_file_arg, TextIOWrapper)
  assert output_file_arg.mode in ["r+", "w"]
  assert output_file_arg.name == output_file


def test_main_input_not_exists(testdir, mock_testclass):
  input_dir = "alphafold"
  output_file = "output.txt"
  ListFiles.list_files = MagicMock()
  with pytest.raises(NotADirectoryError):
    ListFiles.main(["-o", output_file, input_dir])
  ListFiles.list_files.assert_not_called()


def test_list_files_alphafold3(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4", "RPB5_RPB6", "RPB7_RPB8"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold3_files(output, output) for output in alphafold_outputs]
  [shutil.copy(Path(__file__).parent.joinpath("ranking_scores.csv"),
               f"{output}/ranking_scores.csv") for output in alphafold_outputs]
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  for output in alphafold_outputs:
    assert f"{output}/{output}_confidences.json" in files
    assert f"{output}/{output}_data.json" in files
    assert f"{output}/{output}_model.cif" in files
    assert f"{output}/{output}_summary_confidences.json" in files
    assert f"{output}/ranking_scores.csv" in files
    assert f"{output}/TERMS_OF_USE.md" in files


def test_list_files_alphafold(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4", "RPB5_RPB6", "RPB7_RPB8"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB5_RPB6/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB7_RPB8/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/result_model_1_multimer_v3_pred_4.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  assert "RPB3_RPB4/ranked_0.pdb" in files
  assert "RPB3_RPB4/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/result_model_4_multimer_v3_pred_0.pkl" in files
  assert "RPB5_RPB6/ranking_debug.json" in files
  assert "RPB5_RPB6/relax_metrics.json" in files
  assert "RPB5_RPB6/timings.json" in files
  assert "RPB5_RPB6/ranked_0.pdb" in files
  assert "RPB5_RPB6/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB5_RPB6/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB5_RPB6/result_model_1_multimer_v3_pred_4.pkl" in files
  assert "RPB7_RPB8/ranking_debug.json" in files
  assert "RPB7_RPB8/relax_metrics.json" in files
  assert "RPB7_RPB8/timings.json" in files
  assert "RPB7_RPB8/ranked_0.pdb" in files
  assert "RPB7_RPB8/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB7_RPB8/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB7_RPB8/result_model_4_multimer_v3_pred_0.pkl" in files
  assert len(files) == 28


def test_list_files_alphafold_all_pdb(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pdb=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  for j in range(0, 25):
    assert f"RPB1_RPB2/ranked_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB1_RPB2/relaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
      assert f"RPB1_RPB2/unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
  assert "RPB1_RPB2/result_model_1_multimer_v3_pred_4.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  for j in range(0, 25):
    assert f"RPB1_RPB2/ranked_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB3_RPB4/relaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
      assert f"RPB3_RPB4/unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
  assert "RPB3_RPB4/result_model_4_multimer_v3_pred_0.pkl" in files
  assert len(files) == 158


def test_list_files_alphafold_all_pkl(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pkl=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB1_RPB2/result_model_{i}_multimer_v3_pred_{j}.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  assert "RPB3_RPB4/ranked_0.pdb" in files
  assert "RPB3_RPB4/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB3_RPB4/result_model_{i}_multimer_v3_pred_{j}.pkl" in files
  assert len(files) == 62


def test_list_files_alphafold_no_best_pkl(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, best_pkl=False)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB1_RPB2/result_model_{i}_multimer_v3_pred_{j}.pkl" not in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  assert "RPB3_RPB4/ranked_0.pdb" in files
  assert "RPB3_RPB4/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB3_RPB4/result_model_{i}_multimer_v3_pred_{j}.pkl" not in files
  assert len(files) == 12


def test_list_files_alphafold_all_pdb_all_pkl(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pdb=True, all_pkl=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  for j in range(0, 25):
    assert f"RPB1_RPB2/ranked_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB1_RPB2/relaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
      assert f"RPB1_RPB2/unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB1_RPB2/result_model_{i}_multimer_v3_pred_{j}.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  for j in range(0, 25):
    assert f"RPB1_RPB2/ranked_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB3_RPB4/relaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
      assert f"RPB3_RPB4/unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
  for i in range(1, 6):
    for j in range(0, 5):
      assert f"RPB3_RPB4/result_model_{i}_multimer_v3_pred_{j}.pkl" in files
  assert len(files) == 206


def test_list_files_alphafold_missing_ranking(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  os.remove("RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/result_model_1_multimer_v3_pred_4.pkl" in files
  assert len(files) == 7


def test_list_files_af2complex(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4", "RPB5_RPB6", "RPB7_RPB8"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"],
             ["240603_519798"], ["240525_625637"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 4)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 4)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"),
              f"RPB5_RPB6/ranking_all_240603_519798.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB7_RPB8/ranking_all_240525_625637.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/relaxed_model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pkl" in files
  assert f"RPB5_RPB6/ranking_all_{run_ids[2][0]}.json" in files
  assert f"RPB5_RPB6/relax_metrics.json" in files
  assert f"RPB5_RPB6/model_5_multimer_v3_p1_{run_ids[2][0]}.pdb" in files
  assert f"RPB5_RPB6/relaxed_model_5_multimer_v3_p1_{run_ids[2][0]}.pdb" in files
  assert f"RPB5_RPB6/model_5_multimer_v3_p1_{run_ids[2][0]}.pkl" in files
  assert f"RPB7_RPB8/ranking_all_{run_ids[3][0]}.json" in files
  assert f"RPB7_RPB8/relax_metrics.json" in files
  assert f"RPB7_RPB8/model_5_multimer_v3_p1_{run_ids[3][0]}.pdb" in files
  assert f"RPB7_RPB8/relaxed_model_5_multimer_v3_p1_{run_ids[3][0]}.pdb" in files
  assert f"RPB7_RPB8/model_5_multimer_v3_p1_{run_ids[3][0]}.pkl" in files
  assert len(files) == 24


def test_list_files_af2complex_metric(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, metric="pitms")
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  assert f"RPB1_RPB2/model_4_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/relaxed_model_4_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/model_4_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_1_multimer_v3_p1_{run_ids[1][0]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_1_multimer_v3_p1_{run_ids[1][0]}.pdb" in files
  assert f"RPB3_RPB4/model_1_multimer_v3_p1_{run_ids[1][0]}.pkl" in files
  assert len(files) == 14


def test_list_files_af2complex_all_pdb(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pdb=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  for i in range(1, 6):
    assert f"RPB1_RPB2/model_{i}_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
    assert f"RPB1_RPB2/relaxed_model_{i}_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pdb" in files
    assert f"RPB3_RPB4/relaxed_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pdb" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pkl" in files
  assert len(files) == 30


def test_list_files_af2complex_no_best_pkl(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, best_pkl=False)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/relaxed_model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  for i in range(1, 6):
    assert f"RPB1_RPB2/model_{i}_multimer_v3_p1_{run_ids[0][0]}.pkl" not in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pkl" not in files
  assert len(files) == 12


def test_list_files_af2complex_all_pkl(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pkl=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/relaxed_model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  for i in range(1, 6):
    assert f"RPB1_RPB2/model_{i}_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pkl" in files
  assert len(files) == 22


def test_list_files_af2complex_all_pdb_all_pkl(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out, all_pdb=True, all_pkl=True)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  for i in range(1, 6):
    assert f"RPB1_RPB2/model_{i}_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
    assert f"RPB1_RPB2/relaxed_model_{i}_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
    assert f"RPB1_RPB2/model_{i}_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pdb" in files
    assert f"RPB3_RPB4/relaxed_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pdb" in files
    assert f"RPB3_RPB4/model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.pkl" in files
  assert len(files) == 38


def test_list_files_af2complex_missing_ranking(testdir, mock_testclass):
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  os.remove("RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  os.remove("RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  for i in range(1, 6):
    if i == 4:
      assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" not in files
    else:
      assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_3_multimer_v3_p1_{run_ids[1][2]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_3_multimer_v3_p1_{run_ids[1][2]}.pdb" in files
  assert f"RPB3_RPB4/model_3_multimer_v3_p1_{run_ids[1][2]}.pkl" in files
  assert len(files) == 8


def test_list_files_progress(testdir, mock_testclass):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out, patch("tqdm.tqdm",
                                                   return_value=alphafold_outputs) as mock_tqdm:
    ListFiles.list_files("", output_out, progress=True)
    mock_tqdm.assert_called_once_with(ANY)
    assert len(mock_tqdm.call_args.args[0]) == len(alphafold_outputs)
    for alphafold_output in alphafold_outputs:
      assert alphafold_output in mock_tqdm.call_args.args[0]
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/result_model_1_multimer_v3_pred_4.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  assert "RPB3_RPB4/ranked_0.pdb" in files
  assert "RPB3_RPB4/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/result_model_4_multimer_v3_pred_0.pkl" in files
  assert len(files) == 14


def test_list_files_std_out(testdir, mock_testclass, capsys):
  alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB1_RPB2/ranking_debug.json")
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug_2.json"),
              "RPB3_RPB4/ranking_debug.json")
  ListFiles.list_files("", sys.stdout)
  out, err = capsys.readouterr()
  files = out.rstrip('\r\n').split("\n")
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert "RPB1_RPB2/ranking_debug.json" in files
  assert "RPB1_RPB2/relax_metrics.json" in files
  assert "RPB1_RPB2/timings.json" in files
  assert "RPB1_RPB2/ranked_0.pdb" in files
  assert "RPB1_RPB2/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB1_RPB2/result_model_1_multimer_v3_pred_4.pkl" in files
  assert "RPB3_RPB4/ranking_debug.json" in files
  assert "RPB3_RPB4/relax_metrics.json" in files
  assert "RPB3_RPB4/timings.json" in files
  assert "RPB3_RPB4/ranked_0.pdb" in files
  assert "RPB3_RPB4/relaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/unrelaxed_model_4_multimer_v3_pred_0.pdb" in files
  assert "RPB3_RPB4/result_model_4_multimer_v3_pred_0.pkl" in files
  assert len(files) == 14


def test_list_files_alphafold_af2complex_mix(testdir, mock_testclass):
  alphafold_outputs = ["RPB5_RPB6"]
  [testdir.mkdir(output) for output in alphafold_outputs]
  [create_alphafold_files(output) for output in alphafold_outputs]
  shutil.copy(Path(__file__).parent.joinpath("ranking_debug.json"),
              "RPB5_RPB6/ranking_debug.json")
  af2complex_outputs = ["RPB1_RPB2", "RPB3_RPB4"]
  run_ids = [["240525_625637"],
             ["240624_450814", "240624_525588", "240624_260995",
              "240624_214651", "240624_536961"]]
  [testdir.mkdir(af2complex_outputs[i]) for i in range(0, 2)]
  [create_af2complex_files(af2complex_outputs[i], run_ids[i]) for i in
   range(0, 2)]
  shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"),
              f"RPB1_RPB2/ranking_all_240525_625637.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_1_multimer_v3_p1_240624_450814.json"),
      f"RPB3_RPB4/ranking_model_1_multimer_v3_p1_240624_450814.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_2_multimer_v3_p1_240624_525588.json"),
      f"RPB3_RPB4/ranking_model_2_multimer_v3_p1_240624_525588.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_3_multimer_v3_p1_240624_260995.json"),
      f"RPB3_RPB4/ranking_model_3_multimer_v3_p1_240624_260995.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_4_multimer_v3_p1_240624_214651.json"),
      f"RPB3_RPB4/ranking_model_4_multimer_v3_p1_240624_214651.json")
  shutil.copy(Path(__file__).parent.joinpath(
      "ranking_model_5_multimer_v3_p1_240624_536961.json"),
      f"RPB3_RPB4/ranking_model_5_multimer_v3_p1_240624_536961.json")
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    ListFiles.list_files("", output_out)
  with open(output_file, 'r') as output_in:
    files = output_in.readlines()
  files = [file.rstrip('\r\n') for file in files]
  [print(file) for file in files]
  assert f"RPB1_RPB2/ranking_all_{run_ids[0][0]}.json" in files
  assert f"RPB1_RPB2/relax_metrics.json" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/relaxed_model_5_multimer_v3_p1_{run_ids[0][0]}.pdb" in files
  assert f"RPB1_RPB2/model_5_multimer_v3_p1_{run_ids[0][0]}.pkl" in files
  for i in range(1, 6):
    assert f"RPB3_RPB4/ranking_model_{i}_multimer_v3_p1_{run_ids[1][i - 1]}.json" in files
  assert f"RPB3_RPB4/relax_metrics.json" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/relaxed_model_4_multimer_v3_p1_{run_ids[1][3]}.pdb" in files
  assert f"RPB3_RPB4/model_4_multimer_v3_p1_{run_ids[1][3]}.pkl" in files
  assert "RPB5_RPB6/ranking_debug.json" in files
  assert "RPB5_RPB6/relax_metrics.json" in files
  assert "RPB5_RPB6/timings.json" in files
  assert "RPB5_RPB6/ranked_0.pdb" in files
  assert "RPB5_RPB6/relaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB5_RPB6/unrelaxed_model_1_multimer_v3_pred_4.pdb" in files
  assert "RPB5_RPB6/result_model_1_multimer_v3_pred_4.pkl" in files
  assert len(files) == 21


def test_get_best_model_alphafold(testdir, mock_testclass):
  ranking_file = Path(__file__).parent.joinpath("ranking_debug.json")
  best_model = ListFiles.get_best_model([ranking_file])
  assert best_model == "model_1_multimer_v3_pred_4"


def test_get_best_model_af2complex(testdir, mock_testclass):
  ranking_file = Path(__file__).parent.joinpath(
      "ranking_all_240525_625637.json")
  best_model = ListFiles.get_best_model([ranking_file])
  assert best_model == "model_5_multimer_v3_p1_240525_625637"


def test_get_best_model_af2complex_metric(testdir, mock_testclass):
  ranking_file = Path(__file__).parent.joinpath(
      "ranking_all_240525_625637.json")
  best_model = ListFiles.get_best_model([ranking_file], "pitms")
  assert best_model == "model_4_multimer_v3_p1_240525_625637"


def test_get_best_model_af2complex_per_model_ranking(testdir, mock_testclass):
  ranking_files = [Path(__file__).parent.joinpath(
      f"ranking_model_1_multimer_v3_p1_240624_450814.json"),
    Path(__file__).parent.joinpath(
        f"ranking_model_2_multimer_v3_p1_240624_525588.json"),
    Path(__file__).parent.joinpath(
        f"ranking_model_3_multimer_v3_p1_240624_260995.json"),
    Path(__file__).parent.joinpath(
        f"ranking_model_4_multimer_v3_p1_240624_214651.json"),
    Path(__file__).parent.joinpath(
        f"ranking_model_5_multimer_v3_p1_240624_536961.json")]
  best_model = ListFiles.get_best_model(ranking_files)
  assert best_model == "model_4_multimer_v3_p1_240624_214651"


def test_get_best_model_af2complex_recycled_same_score(testdir, mock_testclass):
  ranking_file = Path(__file__).parent.joinpath(
      "ranking_all_240909_072299.json")
  best_model = ListFiles.get_best_model([ranking_file])
  assert best_model == "model_3_multimer_v3_p1_240909_072299"
