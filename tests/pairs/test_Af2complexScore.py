import argparse
import shutil
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import pytest

from pairs import Af2complexScore


@pytest.fixture
def mock_testclass():
    _multi_interaction_score = Af2complexScore.multi_interaction_score
    _parse_rankings = Af2complexScore.parse_rankings
    _parse_mapping = Af2complexScore.parse_mapping
    yield
    Af2complexScore.multi_interaction_score = _multi_interaction_score
    Af2complexScore.parse_rankings = _parse_rankings
    Af2complexScore.parse_mapping = _parse_mapping


def test_main(testdir, mock_testclass):
    Af2complexScore.multi_interaction_score = MagicMock()
    Af2complexScore.main([])
    Af2complexScore.multi_interaction_score.assert_called_once_with(
        input_dir="", output_file=ANY, name=r"([\w-]+)__([\w-]+)",
        metrics=["interface"], progress=False, recycled=False,
        mapping_file=None, source_column=0, converted_column=1)
    output_file = Af2complexScore.multi_interaction_score.call_args.kwargs["output_file"]
    assert isinstance(output_file, TextIOWrapper)
    assert output_file.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
    output = "output.txt"
    metrics = ["pitms", "plddts"]
    name = r"(\w+)_(\w+)"
    mapping = "mapping.txt"
    Path(mapping).touch()
    source_column = 2
    converted_column = 3
    Af2complexScore.multi_interaction_score = MagicMock()
    Af2complexScore.main(["-i", str(testdir), "-o", output, "-m", metrics[0], metrics[1], "-n", name, "-p", "-R",
                          "-M", mapping, "-S", str(source_column + 1), "-C", str(converted_column + 1)])
    Af2complexScore.multi_interaction_score.assert_called_once_with(
        input_dir=str(testdir), output_file=ANY, name=name,
        metrics=metrics, progress=True, recycled=True,
        mapping_file=ANY, source_column=2, converted_column=3)
    output_file = Af2complexScore.multi_interaction_score.call_args.kwargs["output_file"]
    mapping_file = Af2complexScore.multi_interaction_score.call_args.kwargs["mapping_file"]
    assert isinstance(output_file, TextIOWrapper)
    assert output_file.name == output
    assert output_file.mode == "w"
    assert isinstance(mapping_file, TextIOWrapper)
    assert mapping_file.name == mapping
    assert mapping_file.mode == "r"


def test_main_long_parameters(testdir, mock_testclass):
    output = "output.txt"
    metrics = ["pitms", "plddts"]
    name = r"(\w+)_(\w+)"
    mapping = "mapping.txt"
    Path(mapping).touch()
    source_column = 2
    converted_column = 3
    Af2complexScore.multi_interaction_score = MagicMock()
    Af2complexScore.main(["--input", str(testdir), "--output", output, "--metric", metrics[0], metrics[1],
                          "--name", name, "--progress", "--recycled",
                          "--mapping", mapping, "--source_column", str(source_column + 1),
                          "--converted_column", str(converted_column + 1)])
    Af2complexScore.multi_interaction_score.assert_called_once_with(
        input_dir=str(testdir), output_file=ANY, name=name,
        metrics=metrics, progress=True, recycled=True,
        mapping_file=ANY, source_column=2, converted_column=3)
    output_file = Af2complexScore.multi_interaction_score.call_args.kwargs["output_file"]
    mapping_file = Af2complexScore.multi_interaction_score.call_args.kwargs["mapping_file"]
    assert isinstance(output_file, TextIOWrapper)
    assert output_file.name == output
    assert output_file.mode == "w"
    assert isinstance(mapping_file, TextIOWrapper)
    assert mapping_file.name == mapping
    assert mapping_file.mode == "r"


def test_main_no_metrics(testdir, mock_testclass):
    Af2complexScore.multi_interaction_score = MagicMock()
    with pytest.raises(SystemExit):
        Af2complexScore.main(["-m"])
    Af2complexScore.multi_interaction_score.assert_not_called()


def test_multi_interaction_score(testdir, mock_testclass):
    ranking_file_1 = "POLR2A__POLR2B/ranking_all_240525_625637.json"
    ranking_file_2 = "POLR2A__POLR2C/ranking_all_240603_519798.json"
    Path(ranking_file_1).parent.mkdir()
    Path(ranking_file_2).parent.mkdir()
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"), ranking_file_1)
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"), ranking_file_2)
    ranking_1 = Af2complexScore.Ranking("model_5_multimer_v3_p1", 0.7772, 180, 210, 0.7059, 89.2, 0.8952)
    ranking_2 = Af2complexScore.Ranking("model_4_multimer_v3_p1", 0.7601, 86, 93, 0.783, 90.04, 0.8985)
    output = "output.txt"
    Af2complexScore.parse_rankings = MagicMock(side_effect=[ranking_1, ranking_2])
    Af2complexScore.parse_mapping = MagicMock()
    with open(output, "w") as output_out:
        Af2complexScore.multi_interaction_score(output_file=output_out)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_1, "interface", False)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_2, "interface", False)
    Af2complexScore.parse_mapping.assert_not_called()
    with open(output, "r") as output_in:
        assert output_in.readline() == "Bait\tTarget\tInterface score\n"
        assert output_in.readline() == "POLR2A\tPOLR2B\t0.7772\n"
        assert output_in.readline() == "POLR2A\tPOLR2C\t0.7601\n"


def test_multi_interaction_score_parameters(testdir, mock_testclass):
    testdir.mkdir("rankings")
    ranking_file_1 = "rankings/RPB-1___RPB-2/ranking_all_240525_625637.json"
    ranking_file_2 = "rankings/RPB-1___RPB-3/ranking_all_240603_519798.json"
    Path(ranking_file_1).parent.mkdir()
    Path(ranking_file_2).parent.mkdir()
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"), ranking_file_1)
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"), ranking_file_2)
    ranking_1 = Af2complexScore.Ranking("model_5_multimer_v3_p1", 0.7772, 180, 210, 0.7059, 89.2, 0.8952)
    ranking_2 = Af2complexScore.Ranking("model_4_multimer_v3_p1", 0.7601, 86, 93, 0.783, 90.04, 0.8985)
    output = "output.txt"
    metrics = ["ptms", "interface", "plddts", "pitms"]
    mappings_file = "mappings.txt"
    Path(mappings_file).touch()
    mappings = {"RPB-1": "POLR2A", "RPB-2": "POLR2B", "RPB-3": "POLR2C"}
    Af2complexScore.parse_rankings = MagicMock(side_effect=[ranking_1, ranking_2])
    Af2complexScore.parse_mapping = MagicMock(return_value=mappings)
    with open(output, "w") as output_out, open(mappings_file, "r") as mappings_in:
        Af2complexScore.multi_interaction_score("rankings", output_out, r"([\w-]+)___([\w-]+)",
                                                metrics, False, True, mappings_in, 2, 3)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_1, "ptms", True)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_2, "ptms", True)
    Af2complexScore.parse_mapping.assert_called_once_with(ANY, 2, 3)
    mappings_file_arg = Af2complexScore.parse_mapping.call_args.args[0]
    assert mappings_file_arg.name == mappings_file
    assert mappings_file_arg.mode == "r"
    with open(output, "r") as output_in:
        assert output_in.readline() == "Bait\tTarget\tpTM\tInterface score\tpLDDT\tTM score (piTM)\n"
        assert output_in.readline() == "POLR2A\tPOLR2B\t0.8952\t0.7772\t89.2\t0.7059\n"
        assert output_in.readline() == "POLR2A\tPOLR2C\t0.8985\t0.7601\t90.04\t0.783\n"


def test_multi_interaction_score_progress(testdir, mock_testclass):
    ranking_file_1 = "POLR2A__POLR2B/ranking_all_240525_625637.json"
    ranking_file_2 = "POLR2A__POLR2C/ranking_all_240603_519798.json"
    Path(ranking_file_1).parent.mkdir()
    Path(ranking_file_2).parent.mkdir()
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"), ranking_file_1)
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"), ranking_file_2)
    ranking_files = [ranking_file_1, ranking_file_2]
    ranking_1 = Af2complexScore.Ranking("model_5_multimer_v3_p1", 0.7772, 180, 210, 0.7059, 89.2, 0.8952)
    ranking_2 = Af2complexScore.Ranking("model_4_multimer_v3_p1", 0.7601, 86, 93, 0.783, 90.04, 0.8985)
    output = "output.txt"
    Af2complexScore.parse_rankings = MagicMock(side_effect=[ranking_1, ranking_2])
    Af2complexScore.parse_mapping = MagicMock()
    with open(output, "w") as output_out, patch("tqdm.tqdm", return_value=ranking_files) as mock_tqdm:
        Af2complexScore.multi_interaction_score(output_file=output_out, progress=True)
        mock_tqdm.assert_called_once_with(ranking_files)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_1, "interface", False)
    Af2complexScore.parse_rankings.assert_any_call(ranking_file_2, "interface", False)
    Af2complexScore.parse_mapping.assert_not_called()
    with open(output, "r") as output_in:
        assert output_in.readline() == "Bait\tTarget\tInterface score\n"
        assert output_in.readline() == "POLR2A\tPOLR2B\t0.7772\n"
        assert output_in.readline() == "POLR2A\tPOLR2C\t0.7601\n"


def test_multi_interaction_score_empty_metrics(testdir, mock_testclass):
    ranking_file_1 = "POLR2A__POLR2B/ranking_all_240525_625637.json"
    ranking_file_2 = "POLR2A__POLR2C/ranking_all_240603_519798.json"
    Path(ranking_file_1).parent.mkdir()
    Path(ranking_file_2).parent.mkdir()
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"), ranking_file_1)
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"), ranking_file_2)
    ranking_1 = Af2complexScore.Ranking("model_5_multimer_v3_p1", 0.7772, 180, 210, 0.7059, 89.2, 0.8952)
    ranking_2 = Af2complexScore.Ranking("model_4_multimer_v3_p1", 0.7601, 86, 93, 0.783, 90.04, 0.8985)
    output = "output.txt"
    Af2complexScore.parse_rankings = MagicMock(side_effect=[ranking_1, ranking_2])
    Af2complexScore.parse_mapping = MagicMock()
    with pytest.raises(AssertionError):
        with open(output, "w") as output_out:
            Af2complexScore.multi_interaction_score(output_file=output_out, metrics=[])
    Af2complexScore.parse_rankings.assert_not_called()
    Af2complexScore.parse_mapping.assert_not_called()


def test_multi_interaction_score_invalid_metrics(testdir, mock_testclass):
    ranking_file_1 = "POLR2A__POLR2B/ranking_all_240525_625637.json"
    ranking_file_2 = "POLR2A__POLR2C/ranking_all_240603_519798.json"
    Path(ranking_file_1).parent.mkdir()
    Path(ranking_file_2).parent.mkdir()
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240525_625637.json"), ranking_file_1)
    shutil.copy(Path(__file__).parent.joinpath("ranking_all_240603_519798.json"), ranking_file_2)
    ranking_1 = Af2complexScore.Ranking("model_5_multimer_v3_p1", 0.7772, 180, 210, 0.7059, 89.2, 0.8952)
    ranking_2 = Af2complexScore.Ranking("model_4_multimer_v3_p1", 0.7601, 86, 93, 0.783, 90.04, 0.8985)
    output = "output.txt"
    Af2complexScore.parse_rankings = MagicMock(side_effect=[ranking_1, ranking_2])
    Af2complexScore.parse_mapping = MagicMock()
    with pytest.raises(AssertionError):
        with open(output, "w") as output_out:
            Af2complexScore.multi_interaction_score(output_file=output_out, metrics=["test"])
    Af2complexScore.parse_rankings.assert_not_called()
    Af2complexScore.parse_mapping.assert_not_called()


def test_parse_rankings(testdir, mock_testclass):
    ranking_file = Path(__file__).parent.joinpath("ranking_all_240525_625637.json")
    ranking = Af2complexScore.parse_rankings(ranking_file)
    assert ranking.model == "model_5_multimer_v3_p1_240525_625637"
    assert ranking.interaction_score == 0.7772
    assert ranking.interfacial_residue_number == 180
    assert ranking.interfacial_contact_number == 210
    assert ranking.pitms == 0.7059
    assert ranking.plddts == 89.2
    assert ranking.ptms == 0.8952


def test_parse_rankings_metric(testdir, mock_testclass):
    ranking_file = Path(__file__).parent.joinpath("ranking_all_240525_625637.json")
    ranking = Af2complexScore.parse_rankings(ranking_file, "pitms")
    assert ranking.model == "model_4_multimer_v3_p1_240525_625637"
    assert ranking.interaction_score == 0.7655
    assert ranking.interfacial_residue_number == 181
    assert ranking.interfacial_contact_number == 208
    assert ranking.pitms == 0.7073
    assert ranking.plddts == 88.25
    assert ranking.ptms == 0.8937


def test_parse_rankings_metric_use_recycle(testdir, mock_testclass):
    ranking_file = Path(__file__).parent.joinpath("ranking_all_240525_625637.json")
    ranking = Af2complexScore.parse_rankings(ranking_file, "interface", True)
    assert ranking.model == "model_5_multimer_v3_p1_240525_625637_recycled_02"
    assert ranking.interaction_score == 0.7782
    assert ranking.interfacial_residue_number == 181
    assert ranking.interfacial_contact_number == 211
    assert ranking.pitms == 0.7075
    assert ranking.plddts == 89.31
    assert ranking.ptms == 0.8964


def test_parse_mapping(testdir, mock_testclass):
    mapping_file = "mapping_file.txt"
    with open(mapping_file, 'w') as mapping_out:
        mapping_out.write("RPB1_HUMAN\tPOLR2A\n")
        mapping_out.write("NOGENE_HUMAN\t\n")
        mapping_out.write("RPB2_HUMAN\tPOLR2B\n")
    with open(mapping_file, 'r') as mapping_in:
        mappings = Af2complexScore.parse_mapping(mapping_file=mapping_in)
    assert "RPB1_HUMAN" in mappings
    assert mappings["RPB1_HUMAN"] == "POLR2A"
    assert "RPB2_HUMAN" in mappings
    assert mappings["RPB2_HUMAN"] == "POLR2B"
    assert "NOGENE_HUMAN" not in mappings
