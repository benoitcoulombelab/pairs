import os
import shutil
from io import TextIOWrapper
from unittest.mock import MagicMock, ANY, patch
from pathlib import Path
import sys

import pytest

from pairs import ListFiles


@pytest.fixture
def mock_testclass():
    _list_files = ListFiles.list_files
    _relaxation_models = ListFiles.relaxation_models
    _alphafold_archive_files = ListFiles.alphafold_archive_files
    yield
    ListFiles.list_files = _list_files
    ListFiles.relaxation_models = _relaxation_models
    ListFiles.alphafold_archive_files = _alphafold_archive_files


def create_alphafold_files(alphafold_output):
    os.mkdir(os.path.join(alphafold_output, "msas"))
    os.mkdir(os.path.join(alphafold_output, "msas/A"))
    os.mkdir(os.path.join(alphafold_output, "msas/B"))
    create_files = ["features.pkl", "ranking_debug.json", "relax_metrics.json", "timings.json",
                    "msas/chain_id_map.json"]
    create_files.extend([f"ranked_{i}.pdb" for i in range(0, 25)])
    create_files.extend([f"unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" for i in range(1, 6) for j in range(0, 5)])
    create_files.extend([f"relaxed_model_{i}_multimer_v3_pred_{j}.pdb" for i in range(1, 6) for j in range(0, 5)])
    create_files.extend([f"result_model_{i}_multimer_v3_pred_{j}.pkl" for i in range(1, 6) for j in range(0, 5)])
    create_msas_files = ["mgnify_hits.sto", "small_bfd_hits.sto", "uniref90_hits.sto",
                         "pdb_hits.sto", "uniprot_hits.sto"]
    [open(os.path.join(alphafold_output, file), 'w') for file in create_files]
    [open(os.path.join(alphafold_output, "msas/A", file), 'w') for file in create_msas_files]
    [open(os.path.join(alphafold_output, "msas/B", file), 'w') for file in create_msas_files]


def test_main(testdir, mock_testclass):
    ListFiles.list_files = MagicMock()
    ListFiles.main([])
    ListFiles.list_files.assert_called_once_with(
        input_dir="", output_file=ANY, progress=False)
    output_file = ListFiles.list_files.call_args.kwargs["output_file"]
    assert isinstance(output_file, TextIOWrapper)
    assert output_file.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
    input_dir = "alphafold"
    output_file = "output.txt"
    testdir.mkdir(input_dir)
    ListFiles.list_files = MagicMock()
    ListFiles.main(["-o", output_file, "-p", input_dir])
    ListFiles.list_files.assert_called_once_with(
        input_dir=input_dir, output_file=ANY, progress=True)
    output_file_arg = ListFiles.list_files.call_args.kwargs["output_file"]
    assert isinstance(output_file_arg, TextIOWrapper)
    assert output_file_arg.mode in ["r+", "w"]
    assert output_file_arg.name == output_file


def test_main_long_parameters(testdir, mock_testclass):
    input_dir = "alphafold"
    output_file = "output.txt"
    testdir.mkdir(input_dir)
    ListFiles.list_files = MagicMock()
    ListFiles.main(["--output", output_file, "--progress", input_dir])
    ListFiles.list_files.assert_called_once_with(
        input_dir=input_dir, output_file=ANY, progress=True)
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


def test_list_files(testdir, mock_testclass):
    alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4", "RPB5_RPB6", "RPB7_RPB8"]
    alphafold_archive_files_ret = [[f"{output}/ranked_0.pdb", f"{output}/relax_metrics.json"]
                                   for output in alphafold_outputs]
    [testdir.mkdir(output) for output in alphafold_outputs]
    [open(os.path.join(output, "ranked_0.pdb"), 'w') for output in alphafold_outputs]
    output_file = "output.txt"
    ListFiles.alphafold_archive_files = MagicMock(side_effect=alphafold_archive_files_ret)
    with open(output_file, 'w') as output_out:
        ListFiles.list_files("", output_file=output_out)
    for output in alphafold_outputs:
        ListFiles.alphafold_archive_files.assert_any_call(output)
    with open(output_file, 'r') as output_in:
        for files in alphafold_archive_files_ret:
            for file in files:
                assert output_in.readline() == file + "\n"


def test_list_files_progress(testdir, mock_testclass):
    alphafold_outputs = ["RPB1_RPB2", "RPB3_RPB4", "RPB5_RPB6", "RPB7_RPB8"]
    alphafold_archive_files_ret = [[f"{output}/ranked_0.pdb", f"{output}/relax_metrics.json"]
                                   for output in alphafold_outputs]
    [testdir.mkdir(output) for output in alphafold_outputs]
    ranked_0_files = [os.path.join(output, "ranked_0.pdb") for output in alphafold_outputs]
    [open(ranked_0_file, 'w') for ranked_0_file in ranked_0_files]
    output_file = "output.txt"
    ListFiles.alphafold_archive_files = MagicMock(side_effect=alphafold_archive_files_ret)
    with open(output_file, 'w') as output_out, patch("tqdm.tqdm", return_value=ranked_0_files) as mock_tqdm:
        ListFiles.list_files("", output_file=output_out, progress=True)
        mock_tqdm.assert_called_once_with(ANY)
        assert len(mock_tqdm.call_args.args[0]) == len(ranked_0_files)
        for ranked_0_file in ranked_0_files:
            assert ranked_0_file in mock_tqdm.call_args.args[0]
    for output in alphafold_outputs:
        ListFiles.alphafold_archive_files.assert_any_call(output)
    with open(output_file, 'r') as output_in:
        for files in alphafold_archive_files_ret:
            for file in files:
                assert output_in.readline() == file + "\n"


def test_alphafold_archive_files(testdir, mock_testclass):
    create_alphafold_files("")
    relaxation_models = ["model_5_multimer_v3_pred_0"]
    ListFiles.relaxation_models = MagicMock(return_value=relaxation_models)
    files = ListFiles.alphafold_archive_files("")
    ListFiles.relaxation_models.assert_called_once_with("")
    assert len(files) == 31
    assert "ranked_0.pdb" in files
    assert "ranking_debug.json" in files
    assert "relax_metrics.json" in files
    assert "timings.json" in files
    for i in range(1, 6):
        for j in range(0, 4):
            f"unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb" in files
    assert "relaxed_model_5_multimer_v3_pred_0.pdb" in files
    assert "result_model_5_multimer_v3_pred_0.pkl" in files


def test_alphafold_archive_files_missing_files(testdir, mock_testclass, capsys):
    input_dir = "alphafold"
    testdir.mkdir(input_dir)
    create_alphafold_files(input_dir)
    os.remove(os.path.join(input_dir, "ranking_debug.json"))
    os.remove(os.path.join(input_dir, "relax_metrics.json"))
    os.remove(os.path.join(input_dir, "timings.json"))
    os.remove(os.path.join(input_dir, "unrelaxed_model_5_multimer_v3_pred_4.pdb"))
    os.remove(os.path.join(input_dir, "relaxed_model_5_multimer_v3_pred_0.pdb"))
    os.remove(os.path.join(input_dir, "result_model_5_multimer_v3_pred_0.pkl"))
    relaxation_models = ["model_5_multimer_v3_pred_0"]
    ListFiles.relaxation_models = MagicMock(return_value=relaxation_models)
    files = ListFiles.alphafold_archive_files(input_dir)
    ListFiles.relaxation_models.assert_called_once_with(input_dir)
    out, err = capsys.readouterr()
    err_lines = err.split("\n")
    assert err_lines[0] == f"File ranking_debug.json is missing from directory {input_dir}"
    assert err_lines[1] == f"File relax_metrics.json is missing from directory {input_dir}"
    assert err_lines[2] == f"File timings.json is missing from directory {input_dir}"
    assert err_lines[3] == f"File unrelaxed_model_5_multimer_v3_pred_4.pdb is missing from directory {input_dir}"
    assert err_lines[4] == f"File relaxed_model_5_multimer_v3_pred_0.pdb is missing from directory {input_dir}"
    assert err_lines[5] == f"File result_model_5_multimer_v3_pred_0.pkl is missing from directory {input_dir}"
    assert len(files) == 25
    assert os.path.join(input_dir, "ranked_0.pdb") in files
    for i in range(1, 6):
        for j in range(0, 4):
            if i != 5 and j != 4:
                os.path.join(input_dir, f"unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb") in files


def test_relaxation_models_json(testdir, mock_testclass):
    input_dir = "alphafold"
    testdir.mkdir(input_dir)
    relaxation_json = os.path.join(input_dir, "relax_metrics.json")
    shutil.copyfile(Path(__file__).parent.joinpath("relax_metrics.json"), relaxation_json)
    relaxation_models = ListFiles.relaxation_models(input_dir)
    assert len(relaxation_models) == 1
    assert "model_5_multimer_v3_pred_0" in relaxation_models


def test_relaxation_models_json_multiple(testdir, mock_testclass):
    input_dir = "alphafold"
    testdir.mkdir(input_dir)
    relaxation_json = os.path.join(input_dir, "relax_metrics.json")
    shutil.copyfile(Path(__file__).parent.joinpath("relax_metrics_multiple.json"), relaxation_json)
    relaxation_models = ListFiles.relaxation_models(input_dir)
    assert len(relaxation_models) == 2
    assert "model_5_multimer_v3_pred_0" in relaxation_models
    assert "model_1_multimer_v3_pred_4" in relaxation_models


def test_relaxation_models_json_missing(testdir, mock_testclass):
    input_dir = "alphafold"
    testdir.mkdir(input_dir)
    relaxed_files = ["relaxed_model_1_multimer_v3_pred_0.pdb", "relaxed_model_2_multimer_v3_pred_1.pdb"]
    [open(os.path.join(input_dir, f), 'w').close() for f in relaxed_files]
    relaxation_models = ListFiles.relaxation_models(input_dir)
    assert len(relaxation_models) == 2
    assert "model_1_multimer_v3_pred_0" in relaxation_models
    assert "model_2_multimer_v3_pred_1" in relaxation_models
