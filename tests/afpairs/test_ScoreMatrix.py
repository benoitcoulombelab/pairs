import math
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import numpy
import pandas
import pytest

from afpairs import ScoreMatrix


@pytest.fixture
def mock_testclass():
    _score_matrix = ScoreMatrix.score_matrix
    _interaction_matrix = ScoreMatrix.interaction_matrix
    _z_score = ScoreMatrix.z_score
    _unique_matrix = ScoreMatrix.unique_matrix
    _parse_scores = ScoreMatrix.parse_scores
    yield
    ScoreMatrix.score_matrix = _score_matrix
    ScoreMatrix.interaction_matrix = _interaction_matrix
    ScoreMatrix.z_score = _z_score
    ScoreMatrix.unique_matrix = _unique_matrix
    ScoreMatrix.parse_scores = _parse_scores


def test_main(testdir, mock_testclass):
    ScoreMatrix.score_matrix = MagicMock()
    stdin_file = "stdin.txt"
    open(stdin_file, 'w').close()
    with open(stdin_file, 'r') as stdin_in, patch('sys.stdin', stdin_in):
        ScoreMatrix.main()
    ScoreMatrix.score_matrix.assert_called_once_with(
        score_files=ANY, output_file=ANY, zscore=False, unique=False)
    score_files = ScoreMatrix.score_matrix.call_args.kwargs["score_files"]
    assert len(score_files) == 1
    score_in = score_files[0]
    assert isinstance(score_in, TextIOWrapper)
    assert score_in.mode == "r"
    output_out = ScoreMatrix.score_matrix.call_args.kwargs["output_file"]
    assert isinstance(output_out, TextIOWrapper)
    assert output_out.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
    scores_file_1 = "scores1.txt"
    scores_file_2 = "scores2.txt"
    open(scores_file_1, 'w').close()
    open(scores_file_2, 'w').close()
    output_file = "output.txt"
    ScoreMatrix.score_matrix = MagicMock()
    ScoreMatrix.main(["-z", "-u", "-o", output_file, scores_file_1, scores_file_2])
    ScoreMatrix.score_matrix.assert_called_once_with(
        score_files=ANY, output_file=ANY, zscore=True, unique=True)
    score_files = ScoreMatrix.score_matrix.call_args.kwargs["score_files"]
    assert len(score_files) == 2
    score_in = score_files[0]
    assert score_in.name == scores_file_1
    assert score_in.mode == "r"
    score_in = score_files[1]
    assert score_in.name == scores_file_2
    assert score_in.mode == "r"
    output_out = ScoreMatrix.score_matrix.call_args.kwargs["output_file"]
    assert output_out.name == output_file
    assert output_out.mode == "w"


def test_main_long_parameters(testdir, mock_testclass):
    scores_file_1 = "scores1.txt"
    scores_file_2 = "scores2.txt"
    open(scores_file_1, 'w').close()
    open(scores_file_2, 'w').close()
    output_file = "output.txt"
    ScoreMatrix.score_matrix = MagicMock()
    ScoreMatrix.main(["--zscore", "--unique", "--output", output_file, scores_file_1, scores_file_2])
    ScoreMatrix.score_matrix.assert_called_once_with(
        score_files=ANY, output_file=ANY, zscore=True, unique=True)
    score_files = ScoreMatrix.score_matrix.call_args.kwargs["score_files"]
    assert len(score_files) == 2
    score_in = score_files[0]
    assert score_in.name == scores_file_1
    assert score_in.mode == "r"
    score_in = score_files[1]
    assert score_in.name == scores_file_2
    assert score_in.mode == "r"
    output_out = ScoreMatrix.score_matrix.call_args.kwargs["output_file"]
    assert output_out.name == output_file
    assert output_out.mode == "w"


def test_score_matrix(testdir, mock_testclass):
    scores_file_1 = "scores1.txt"
    scores_file_2 = "scores2.txt"
    open(scores_file_1, 'w').close()
    open(scores_file_2, 'w').close()
    output_file = "output.txt"
    interactions_1 = list(range(0, 10))
    interactions_2 = list(range(8, 18))
    ScoreMatrix.parse_scores = MagicMock(side_effect=[interactions_1, interactions_2])
    random_data = numpy.random.rand(5, 3)
    matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C"],
                              index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    ScoreMatrix.interaction_matrix = MagicMock(return_value=matrix)
    print(random_data[0])
    with open(scores_file_1, 'r') as scores_in_1, open(scores_file_2, 'r') as scores_in_2, \
            open(output_file, 'w') as output_out:
        ScoreMatrix.score_matrix([scores_in_1, scores_in_2], output_file=output_out)
        ScoreMatrix.parse_scores.assert_any_call(scores_in_1)
        ScoreMatrix.parse_scores.assert_any_call(scores_in_2)
    ScoreMatrix.interaction_matrix.assert_called_once_with(interactions_1 + interactions_2, zscore=False, unique=False)
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "Target\tPOLR2A\tPOLR2B\tPOLR2C\n"
        assert output_in.readline() == "POLR2A\t" + "\t".join([str(d) for d in random_data[0]]) + "\n"
        assert output_in.readline() == "POLR2B\t" + "\t".join([str(d) for d in random_data[1]]) + "\n"
        assert output_in.readline() == "POLR2C\t" + "\t".join([str(d) for d in random_data[2]]) + "\n"
        assert output_in.readline() == "POLR2D\t" + "\t".join([str(d) for d in random_data[3]]) + "\n"
        assert output_in.readline() == "POLR2E\t" + "\t".join([str(d) for d in random_data[4]]) + "\n"


def test_score_matrix_unique_zscore(testdir, mock_testclass):
    scores_file_1 = "scores1.txt"
    scores_file_2 = "scores2.txt"
    open(scores_file_1, 'w').close()
    open(scores_file_2, 'w').close()
    output_file = "output.txt"
    interactions_1 = list(range(0, 10))
    interactions_2 = list(range(8, 18))
    ScoreMatrix.parse_scores = MagicMock(side_effect=[interactions_1, interactions_2])
    random_data = numpy.random.rand(5, 3)
    matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C"],
                              index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    ScoreMatrix.interaction_matrix = MagicMock(return_value=matrix)
    with open(scores_file_1, 'r') as scores_in_1, open(scores_file_2, 'r') as scores_in_2, \
            open(output_file, 'w') as output_out:
        ScoreMatrix.score_matrix([scores_in_1, scores_in_2], output_file=output_out, zscore=True, unique=True)
        ScoreMatrix.parse_scores.assert_any_call(scores_in_1)
        ScoreMatrix.parse_scores.assert_any_call(scores_in_2)
    ScoreMatrix.interaction_matrix.assert_called_once_with(interactions_1 + interactions_2, zscore=True, unique=True)
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "Target\tPOLR2A\tPOLR2B\tPOLR2C\n"
        assert output_in.readline() == "POLR2A\t" + "\t".join([str(d) for d in random_data[0]]) + "\n"
        assert output_in.readline() == "POLR2B\t" + "\t".join([str(d) for d in random_data[1]]) + "\n"
        assert output_in.readline() == "POLR2C\t" + "\t".join([str(d) for d in random_data[2]]) + "\n"
        assert output_in.readline() == "POLR2D\t" + "\t".join([str(d) for d in random_data[3]]) + "\n"
        assert output_in.readline() == "POLR2E\t" + "\t".join([str(d) for d in random_data[4]]) + "\n"


def test_interaction_matrix(mock_testclass):
    scores_file = Path(__file__).parent.joinpath("polr2_scores.tsv")
    with open(scores_file, 'r') as scores_in:
        interactions = ScoreMatrix.parse_scores(scores_in)
    matrix = ScoreMatrix.interaction_matrix(interactions)
    assert matrix.shape == (12, 12)
    assert list(matrix.columns) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert list(matrix.index) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert math.isnan(matrix.at["POLR2A", "POLR2B"])
    assert abs(matrix.at["POLR2B", "POLR2A"] - 31.48) < 0.01
    assert abs(matrix.at["POLR2B", "POLR2E"] - 0.1736751131457573) < 0.01
    assert abs(matrix.at["POLR2E", "POLR2B"] - 2.141993062131007) < 0.01


def test_interaction_matrix_unique(mock_testclass):
    scores_file = Path(__file__).parent.joinpath("polr2_scores.tsv")
    with open(scores_file, 'r') as scores_in:
        interactions = ScoreMatrix.parse_scores(scores_in)
    random_data = numpy.random.rand(5, 5)
    random_matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                     index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    ScoreMatrix.unique_matrix = MagicMock(return_value=random_matrix)
    matrix = ScoreMatrix.interaction_matrix(interactions, unique=True)
    ScoreMatrix.unique_matrix.assert_called_once_with(ANY)
    before_unique = ScoreMatrix.unique_matrix.call_args.args[0]
    assert before_unique.shape == (12, 12)
    assert list(before_unique.columns) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert list(before_unique.index) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert math.isnan(before_unique.at["POLR2A", "POLR2B"])
    assert abs(before_unique.at["POLR2B", "POLR2A"] - 31.48) < 0.01
    assert abs(before_unique.at["POLR2B", "POLR2E"] - 0.1736751131457573) < 0.01
    assert abs(before_unique.at["POLR2E", "POLR2B"] - 2.141993062131007) < 0.01
    assert matrix.equals(random_matrix)


def test_interaction_matrix_zscore(mock_testclass):
    scores_file = Path(__file__).parent.joinpath("polr2_scores.tsv")
    with open(scores_file, 'r') as scores_in:
        interactions = ScoreMatrix.parse_scores(scores_in)
    random_data = numpy.random.rand(5, 5)
    random_matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                     index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    ScoreMatrix.z_score = MagicMock(return_value=random_matrix)
    matrix = ScoreMatrix.interaction_matrix(interactions, zscore=True)
    ScoreMatrix.z_score.assert_called_once_with(ANY)
    before_zscore = ScoreMatrix.z_score.call_args.args[0]
    assert before_zscore.shape == (12, 12)
    assert list(before_zscore.columns) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert list(before_zscore.index) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert math.isnan(before_zscore.at["POLR2A", "POLR2B"])
    assert abs(before_zscore.at["POLR2B", "POLR2A"] - 31.48) < 0.01
    assert abs(before_zscore.at["POLR2B", "POLR2E"] - 0.1736751131457573) < 0.01
    assert abs(before_zscore.at["POLR2E", "POLR2B"] - 2.141993062131007) < 0.01
    assert matrix.equals(random_matrix)


def test_interaction_matrix_zscore_unique(mock_testclass):
    scores_file = Path(__file__).parent.joinpath("polr2_scores.tsv")
    with open(scores_file, 'r') as scores_in:
        interactions = ScoreMatrix.parse_scores(scores_in)
    random_matrix_1 = pandas.DataFrame(numpy.random.rand(5, 5),
                                       columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                       index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    random_matrix_2 = pandas.DataFrame(numpy.random.rand(5, 5),
                                       columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                       index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    random_matrix_3 = pandas.DataFrame(numpy.random.rand(5, 5),
                                       columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                       index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    ScoreMatrix.unique_matrix = MagicMock(side_effect=[random_matrix_1, random_matrix_2])
    ScoreMatrix.z_score = MagicMock(return_value=random_matrix_3)
    matrix = ScoreMatrix.interaction_matrix(interactions, zscore=True, unique=True)
    assert 2 == ScoreMatrix.unique_matrix.call_count
    unique_matrix_1 = ScoreMatrix.unique_matrix.call_args_list[0].args[0]
    unique_matrix_2 = ScoreMatrix.unique_matrix.call_args_list[1].args[0]
    ScoreMatrix.z_score.assert_called_once_with(ANY)
    zscore_matrix_1 = ScoreMatrix.z_score.call_args.args[0]
    assert unique_matrix_1.shape == (12, 12)
    assert list(unique_matrix_1.columns) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert list(unique_matrix_1.index) == ["POLR2" + c for c in "ABCDEFGHIJKL"]
    assert math.isnan(unique_matrix_1.at["POLR2A", "POLR2B"])
    assert abs(unique_matrix_1.at["POLR2B", "POLR2A"] - 31.48) < 0.01
    assert abs(unique_matrix_1.at["POLR2B", "POLR2E"] - 0.1736751131457573) < 0.01
    assert abs(unique_matrix_1.at["POLR2E", "POLR2B"] - 2.141993062131007) < 0.01
    assert zscore_matrix_1.equals(random_matrix_1)
    assert unique_matrix_2.equals(random_matrix_3)
    assert matrix.equals(random_matrix_2)


def test_zscore(mock_testclass):
    random_matrix = pandas.DataFrame(numpy.random.rand(5, 5),
                                     columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                     index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    means = random_matrix.mean()
    standard_deviations = random_matrix.std(ddof=0)
    matrix = ScoreMatrix.z_score(random_matrix)
    for i in matrix.index:
        for j in matrix.columns:
            assert (random_matrix.at[i, j] - means[j]) / standard_deviations[j] == matrix.at[i, j]


def test_unique_matrix(mock_testclass):
    random_data = numpy.random.rand(5, 5)
    random_matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                     index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    matrix = ScoreMatrix.unique_matrix(random_matrix)
    for i in matrix.index:
        for j in matrix.columns:
            assert max(random_matrix.at[i, j], random_matrix.at[j, i]) == matrix.at[i, j]


def test_unique_matrix_nan(mock_testclass):
    random_data = numpy.random.rand(5, 5)
    random_matrix = pandas.DataFrame(random_data, columns=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
                                     index=["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    random_matrix.loc["POLR2A", "POLR2B"] = math.nan
    random_matrix.loc["POLR2B", "POLR2A"] = math.nan
    matrix = ScoreMatrix.unique_matrix(random_matrix)
    for i in matrix.index:
        for j in matrix.columns:
            if math.isnan(random_matrix.at[i, j]) and math.isnan(random_matrix.at[j, i]):
                assert math.isnan(matrix.at[i, j])
            else:
                assert max(random_matrix.at[i, j], random_matrix.at[j, i]) == matrix.at[i, j]


def test_parse_scores(mock_testclass):
    scores_file = Path(__file__).parent.joinpath("scores.tsv")
    with open(scores_file, 'r') as scores_in:
        interactions = ScoreMatrix.parse_scores(scores_in)
    assert len(interactions) == 5
    interaction = interactions[0]
    assert interaction.bait == "POLR2E"
    assert interaction.target == "POLR2K"
    assert abs(interaction.score - 3.74) < 0.01
    interaction = interactions[1]
    assert interaction.bait == "POLR2E"
    assert interaction.target == "POLR2L"
    assert abs(interaction.score - 3.00) < 0.01
    interaction = interactions[2]
    assert interaction.bait == "POLR2E"
    assert interaction.target == "POLR2A"
    assert abs(interaction.score - 7.54) < 0.01
    interaction = interactions[3]
    assert interaction.bait == "POLR2F"
    assert interaction.target == "POLR2E"
    assert abs(interaction.score - 1.24) < 0.01
    interaction = interactions[4]
    assert interaction.bait == "POLR2H"
    assert interaction.target == "POLR2F"
    assert abs(interaction.score - 4.14) < 0.01
