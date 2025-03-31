import os.path
import shutil
import sys
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import pytest
import smokesignal

from pairs import InteractionScore, MultiInteractionScore


@pytest.fixture
def mock_testclass():
  _file_path = MultiInteractionScore.file_path
  _multi_interaction_score = MultiInteractionScore.multi_interaction_score
  _alphafold_statistics = MultiInteractionScore.alphafold_statistics
  _parse_mapping = MultiInteractionScore.parse_mapping
  _interaction_score = InteractionScore.interaction_score
  _smokesignal_on = smokesignal.on
  yield
  MultiInteractionScore.file_path = _file_path
  MultiInteractionScore.multi_interaction_score = _multi_interaction_score
  MultiInteractionScore.alphafold_statistics = _alphafold_statistics
  MultiInteractionScore.parse_mapping = _parse_mapping
  InteractionScore.interaction_score = _interaction_score
  smokesignal.on = _smokesignal_on


def test_main(testdir, mock_testclass):
  pdbs = ["1.pdb", "2.pdb", "3.pdb"]
  [open(f, 'w').close() for f in pdbs]
  MultiInteractionScore.multi_interaction_score = MagicMock()
  smokesignal.on = MagicMock()
  MultiInteractionScore.main(pdbs)
  MultiInteractionScore.multi_interaction_score.assert_called_once_with(
      input_files=pdbs, output_file=ANY, name=r"([\w-]+)__([\w-]+)",
      stats=False,
      radius=6.0, weight=False, count=False, progress=False,
      first_chains=["A"], second_chains=["B"], partial=False,
      mapping_file=None, source_column=0, converted_column=1)
  output_file = MultiInteractionScore.multi_interaction_score.call_args.kwargs[
    "output_file"]
  assert isinstance(output_file, TextIOWrapper)
  assert output_file.mode in ["r+", "w"]
  smokesignal.on.assert_called_once_with(
      InteractionScore.MISSING_CHAIN_EVENT, InteractionScore.warn_missing_chain,
      max_calls=1)


def test_main_parameters(testdir, mock_testclass):
  pdbs = ["1.pdb", "2.pdb", "3.pdb"]
  [open(f, 'w').close() for f in pdbs]
  mapping_file = "mappings.txt"
  open(mapping_file, 'w').close()
  output_file = "output.txt"
  MultiInteractionScore.multi_interaction_score = MagicMock()
  smokesignal.on = MagicMock()
  MultiInteractionScore.main(
      ["-o", output_file, "-n", r"(\w+)____(\w+)", "-s", "-A", "A,B", "-B",
       "C,D", "-r", "9", "-w", "-c", "-p", "-P",
       "-M", mapping_file, "-S", "2", "-C", "3"] + pdbs)
  MultiInteractionScore.multi_interaction_score.assert_called_once_with(
      input_files=pdbs, output_file=ANY, name=r"(\w+)____(\w+)",
      stats=True,
      radius=9.0, weight=True, count=True, progress=True,
      first_chains=["A", "B"], second_chains=["C", "D"], partial=True,
      mapping_file=ANY, source_column=1, converted_column=2)
  mapping_in = MultiInteractionScore.multi_interaction_score.call_args.kwargs[
    "mapping_file"]
  assert mapping_in.name == mapping_file
  assert mapping_in.mode == "r"
  output_out = MultiInteractionScore.multi_interaction_score.call_args.kwargs[
    "output_file"]
  assert output_out.name == output_file
  assert output_out.mode == "w"
  smokesignal.on.assert_called_once_with(
      InteractionScore.MISSING_CHAIN_EVENT, InteractionScore.warn_missing_chain,
      max_calls=1)


def test_main_long_parameters(testdir, mock_testclass):
  pdbs = ["1.pdb", "2.pdb", "3.pdb"]
  [open(f, 'w').close() for f in pdbs]
  mapping_file = "mappings.txt"
  open(mapping_file, 'w').close()
  output_file = "output.txt"
  MultiInteractionScore.multi_interaction_score = MagicMock()
  smokesignal.on = MagicMock()
  MultiInteractionScore.main(
      ["--output", output_file, "--name", r"(\w+)____(\w+)", "--stats",
       "--first", "A,B", "--second", "C,D",
       "--radius", "9", "--weight", "--count", "--progress", "--partial",
       "--mapping", mapping_file, "--source_column", "2", "--converted_column",
       "3"] + pdbs)
  MultiInteractionScore.multi_interaction_score.assert_called_once_with(
      input_files=pdbs, output_file=ANY, name=r"(\w+)____(\w+)",
      stats=True,
      radius=9.0, weight=True, count=True, progress=True,
      first_chains=["A", "B"], second_chains=["C", "D"], partial=True,
      mapping_file=ANY, source_column=1, converted_column=2)
  mapping_in = MultiInteractionScore.multi_interaction_score.call_args.kwargs[
    "mapping_file"]
  assert mapping_in.name == mapping_file
  assert mapping_in.mode == "r"
  output_out = MultiInteractionScore.multi_interaction_score.call_args.kwargs[
    "output_file"]
  assert output_out.name == output_file
  assert output_out.mode == "w"
  smokesignal.on.assert_called_once_with(
      InteractionScore.MISSING_CHAIN_EVENT, InteractionScore.warn_missing_chain,
      max_calls=1)


def test_multi_interaction_score(testdir, mock_testclass):
  pdbs = ["POLR2A__POLR2B.pdb", "POLR2C__POLR2J-I.pdb", "POLR2D-E__POLR2G.pdb"]
  [open(f, 'w').close() for f in pdbs]
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  with open(output_file, 'w') as output_out:
    MultiInteractionScore.multi_interaction_score(input_files=pdbs,
                                                  output_file=output_out)
  InteractionScore.interaction_score.assert_any_call(pdb=ANY, radius=6.0,
                                                     weight=False, count=False,
                                                     first_chains=["A"],
                                                     second_chains=["B"],
                                                     partial=False)
  for i in range(0, len(pdbs)):
    assert InteractionScore.interaction_score.call_args_list[i].kwargs[
             "pdb"].name == pdbs[i]
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\tScore\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t2.3\n"
    assert output_in.readline() == "POLR2C\tPOLR2J-I\t4.5\n"
    assert output_in.readline() == "POLR2D-E\tPOLR2G\t8.2\n"


def test_multi_interaction_score_parameters(testdir, mock_testclass):
  pdbs = ["POLR2A__POLR2B.pdb", "POLR2C__POLR2J-I.pdb", "POLR2D-E__POLR2G.pdb"]
  [open(f, 'w').close() for f in pdbs]
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  with open(output_file, 'w') as output_out:
    MultiInteractionScore.multi_interaction_score(
        input_files=pdbs, output_file=output_out, radius=9.0, weight=True,
        count=True,
        first_chains=["A", "B"], second_chains=["C", "D"], partial=True)
  InteractionScore.interaction_score.assert_any_call(pdb=ANY, radius=9.0,
                                                     weight=True, count=True,
                                                     first_chains=["A", "B"],
                                                     second_chains=["C", "D"],
                                                     partial=True)
  for i in range(0, len(pdbs)):
    assert InteractionScore.interaction_score.call_args_list[i].kwargs[
             "pdb"].name == pdbs[i]
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\tScore\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t2.3\n"
    assert output_in.readline() == "POLR2C\tPOLR2J-I\t4.5\n"
    assert output_in.readline() == "POLR2D-E\tPOLR2G\t8.2\n"


def test_multi_interaction_score_mapping(testdir, mock_testclass):
  pdbs = ["RPB1_HUMAN__RPB2_HUMAN.pdb", "RPB3_HUMAN__RPB11_HUMAN.pdb",
          "RPB4_HUMAN__RPB7_HUMAN.pdb"]
  [open(f, 'w').close() for f in pdbs]
  mapping_file = "mapping_file.txt"
  open(mapping_file, 'w').close()
  mappings = {"RPB1_HUMAN": "POLR2A", "RPB2_HUMAN": "POLR2B",
              "RPB3_HUMAN": "POLR2C", "RPB7_HUMAN": "POLR2G"}
  MultiInteractionScore.parse_mapping = MagicMock(return_value=mappings)
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  with open(output_file, 'w') as output_out, open(mapping_file,
                                                  'r') as mapping_in:
    MultiInteractionScore.multi_interaction_score(input_files=pdbs,
                                                  output_file=output_out,
                                                  mapping_file=mapping_in)
    MultiInteractionScore.parse_mapping.assert_called_once_with(mapping_in, 0,
                                                                1)
  InteractionScore.interaction_score.assert_any_call(pdb=ANY, radius=6.0,
                                                     weight=False, count=False,
                                                     first_chains=["A"],
                                                     second_chains=["B"],
                                                     partial=False)
  for i in range(0, len(pdbs)):
    assert InteractionScore.interaction_score.call_args_list[i].kwargs[
             "pdb"].name == pdbs[i]
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\tScore\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t2.3\n"
    assert output_in.readline() == "POLR2C\tRPB11_HUMAN\t4.5\n"
    assert output_in.readline() == "RPB4_HUMAN\tPOLR2G\t8.2\n"


def test_multi_interaction_score_invalid_name(testdir, mock_testclass):
  pdbs = ["POLR2A__POLR2B.pdb", "POLR2C__POLR2J.pdb", "POLR2D__POLR2G.pdb"]
  [open(f, 'w').close() for f in pdbs]
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  try:
    with open(output_file, 'w') as output_out:
      MultiInteractionScore.multi_interaction_score(input_files=pdbs,
                                                    name=r"(\w+)",
                                                    output_file=output_out)
  except (AssertionError, IndexError):
    assert True
  else:
    assert False, "Expected AssertionError"


def test_multi_interaction_score_invalid_filename(testdir, mock_testclass):
  pdbs = ["POLR2A.pdb", "POLR2C__POLR2J.pdb", "POLR2D__POLR2G.pdb"]
  [open(f, 'w').close() for f in pdbs]
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  try:
    with open(output_file, 'w') as output_out:
      MultiInteractionScore.multi_interaction_score(input_files=pdbs,
                                                    output_file=output_out)
  except (AssertionError, IndexError):
    assert True
  else:
    assert False, "Expected AssertionError"


def test_multi_interaction_score_unused_chain(testdir, mock_testclass, capsys):
  smokesignal.on(InteractionScore.MISSING_CHAIN_EVENT,
                 InteractionScore.warn_missing_chain, max_calls=1)
  pdb = Path(__file__).parent.joinpath("FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  MultiInteractionScore.multi_interaction_score(input_files=[str(pdb)])
  out, err = capsys.readouterr()
  sys.stderr.write(err)
  assert err == "Chain C present in PDB but not used for scoring\n"


def test_multi_interaction_score_unused_chain_partial(testdir, mock_testclass,
    capsys):
  smokesignal.on(InteractionScore.MISSING_CHAIN_EVENT,
                 InteractionScore.warn_missing_chain, max_calls=1)
  pdb = Path(__file__).parent.joinpath("FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  MultiInteractionScore.multi_interaction_score(input_files=[str(pdb)],
                                                partial=True)
  out, err = capsys.readouterr()
  sys.stderr.write(err)
  assert "Chain C present in PDB but not used for scoring\n" not in err


def test_multi_interaction_score_unused_chain_multiple_pdb(testdir,
    mock_testclass, capsys):
  smokesignal.on(InteractionScore.MISSING_CHAIN_EVENT,
                 InteractionScore.warn_missing_chain, max_calls=1)
  pdb = Path(__file__).parent.joinpath("FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  MultiInteractionScore.multi_interaction_score(
      input_files=[str(pdb), str(pdb)])
  out, err = capsys.readouterr()
  sys.stderr.write(err)
  assert err == "Chain C present in PDB but not used for scoring\n"


def test_multi_interaction_score_progress(testdir, mock_testclass):
  pdbs = ["POLR2A__POLR2B.pdb", "POLR2C__POLR2J-I.pdb", "POLR2D-E__POLR2G.pdb"]
  [open(f, 'w').close() for f in pdbs]
  output_file = "output.txt"
  InteractionScore.interaction_score = MagicMock(side_effect=[2.3, 4.5, 8.2])
  with open(output_file, 'w') as output_out, patch("tqdm.tqdm",
                                                   return_value=pdbs) as mock_tqdm:
    MultiInteractionScore.multi_interaction_score(input_files=pdbs,
                                                  output_file=output_out,
                                                  progress=True)
    mock_tqdm.assert_called_once_with(pdbs)
  InteractionScore.interaction_score.assert_any_call(pdb=ANY, radius=6.0,
                                                     weight=False, count=False,
                                                     first_chains=["A"],
                                                     second_chains=["B"],
                                                     partial=False)
  for i in range(0, len(pdbs)):
    assert InteractionScore.interaction_score.call_args_list[i].kwargs[
             "pdb"].name == pdbs[i]
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\tScore\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t2.3\n"
    assert output_in.readline() == "POLR2C\tPOLR2J-I\t4.5\n"
    assert output_in.readline() == "POLR2D-E\tPOLR2G\t8.2\n"


def test_multi_interaction_score_statistics(testdir, mock_testclass):
  directories = ["POLR2A__POLR2B", "POLR2C__POLR2J", "POLR2D__POLR2G"]
  ranked0_files = [os.path.join(directory, "ranked_0.pdb") for directory in
                   directories]
  [testdir.mkdir(directory) for directory in directories]
  [open(f, 'w').close() for f in ranked0_files]
  as1 = MultiInteractionScore.AlphafoldStatistics(directories[0], 8.2, 0.8)
  as2 = MultiInteractionScore.AlphafoldStatistics(directories[1], 6.5, 0.5)
  as3 = MultiInteractionScore.AlphafoldStatistics(directories[2], 10.3, 0.31)
  MultiInteractionScore.alphafold_statistics = MagicMock(
      side_effect=[as1, as2, as3])
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    MultiInteractionScore.multi_interaction_score(ranked0_files,
                                                  output_file=output_out,
                                                  stats=True)
  for directory in directories:
    MultiInteractionScore.alphafold_statistics.assert_any_call(
        directory=directory, radius=6.0, weight=False, count=False,
        first_chains=["A"], second_chains=["B"], partial=False)
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\t" \
                                   "PAIRS score\tConfidence\tScore * Confidence\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t8.2\t0.8\t6.56\n"
    assert output_in.readline() == "POLR2C\tPOLR2J\t6.5\t0.5\t3.25\n"
    assert output_in.readline() == "POLR2D\tPOLR2G\t10.3\t0.31\t3.193\n"


def test_multi_interaction_score_statistics_parameters(testdir, mock_testclass):
  directories = ["POLR2A__POLR2B", "POLR2C__POLR2J", "POLR2D__POLR2G"]
  ranked0_files = [os.path.join(directory, "ranked_0.pdb") for directory in
                   directories]
  [testdir.mkdir(directory) for directory in directories]
  [open(f, 'w').close() for f in ranked0_files]
  as1 = MultiInteractionScore.AlphafoldStatistics(directories[0], 8.2, 0.8)
  as2 = MultiInteractionScore.AlphafoldStatistics(directories[1], 6.5, 0.5)
  as3 = MultiInteractionScore.AlphafoldStatistics(directories[2], 10.3, 0.31)
  MultiInteractionScore.alphafold_statistics = MagicMock(
      side_effect=[as1, as2, as3])
  output_file = "output.txt"
  with open(output_file, 'w') as output_out:
    MultiInteractionScore.multi_interaction_score(
        ranked0_files, output_file=output_out, stats=True, radius=9.0,
        weight=True, count=True,
        first_chains=["A", "B"], second_chains=["C", "D"], partial=True)
  for directory in directories:
    MultiInteractionScore.alphafold_statistics.assert_any_call(
        directory=directory, radius=9.0, weight=True, count=True,
        first_chains=["A", "B"], second_chains=["C", "D"], partial=True)
  assert os.path.isfile(output_file)
  with open(output_file, 'r') as output_in:
    assert output_in.readline() == "Bait\tTarget\t" \
                                   "PAIRS score\tConfidence\tScore * Confidence\n"
    assert output_in.readline() == "POLR2A\tPOLR2B\t8.2\t0.8\t6.56\n"
    assert output_in.readline() == "POLR2C\tPOLR2J\t6.5\t0.5\t3.25\n"
    assert output_in.readline() == "POLR2D\tPOLR2G\t10.3\t0.31\t3.193\n"


def test_alphafold_statistics(testdir, mock_testclass):
  directory = testdir.mkdir("FAB_5_3__HVM13_MOUSE")
  ranked0_file = Path(__file__).parent.joinpath(
      "FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  shutil.copy(ranked0_file, os.path.join(directory, "ranked_0.pdb"))
  rankings_file = Path(__file__).parent.joinpath("ranking_debug.json")
  shutil.copy(rankings_file, directory)
  unrelaxed_files = ["unrelaxed_model_1_multimer_v3_pred_0.pdb",
                     "unrelaxed_model_1_multimer_v3_pred_1.pdb",
                     "unrelaxed_model_2_multimer_v3_pred_0.pdb"]
  [open(os.path.join(directory, f), 'w').close() for f in unrelaxed_files]
  InteractionScore.interaction_score = MagicMock(
      side_effect=[8.2, 2.3, 4.5, 8.2])
  alphafold_statistics = MultiInteractionScore.alphafold_statistics(directory)
  InteractionScore.interaction_score.assert_any_call(
      pdb=ANY, radius=6.0, weight=False, count=False, first_chains=["A"],
      second_chains=["B"], partial=False)
  assert InteractionScore.interaction_score.call_args_list[0].kwargs["pdb"].name \
         == os.path.join(directory, "ranked_0.pdb")
  assert InteractionScore.interaction_score.call_args_list[0].kwargs[
           "first_chains"] == ["A"]
  assert InteractionScore.interaction_score.call_args_list[0].kwargs[
           "second_chains"] == ["B"]
  assert alphafold_statistics.directory == directory
  assert alphafold_statistics.ranked0_score == 8.2
  assert abs(alphafold_statistics.ranked0_confidence - 0.4154) < 0.0001, \
    f"{alphafold_statistics.ranked0_confidence}"


def test_alphafold_statistics_ranked0_only(testdir, mock_testclass):
  directory = testdir.mkdir("FAB_5_3__HVM13_MOUSE")
  ranked0_file = Path(__file__).parent.joinpath(
      "FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  shutil.copy(ranked0_file, os.path.join(directory, "ranked_0.pdb"))
  InteractionScore.interaction_score = MagicMock(return_value=8.2)
  alphafold_statistics = MultiInteractionScore.alphafold_statistics(directory)
  InteractionScore.interaction_score.assert_called_once_with(
      pdb=ANY, radius=6.0, weight=False, count=False, first_chains=["A"],
      second_chains=["B"], partial=False)
  assert InteractionScore.interaction_score.call_args.kwargs[
           "pdb"].name == os.path.join(directory, "ranked_0.pdb")
  assert alphafold_statistics.directory == directory
  assert alphafold_statistics.ranked0_score == 8.2
  assert alphafold_statistics.ranked0_confidence is None


def test_parse_mapping(testdir, mock_testclass):
  mapping_file = "mapping_file.txt"
  with open(mapping_file, 'w') as mapping_out:
    mapping_out.write("RPB1_HUMAN\tPOLR2A\n")
    mapping_out.write("NOGENE_HUMAN\t\n")
    mapping_out.write("RPB2_HUMAN\tPOLR2B\n")
  with open(mapping_file, 'r') as mapping_in:
    mappings = MultiInteractionScore.parse_mapping(mapping_file=mapping_in)
  assert "RPB1_HUMAN" in mappings
  assert mappings["RPB1_HUMAN"] == "POLR2A"
  assert "RPB2_HUMAN" in mappings
  assert mappings["RPB2_HUMAN"] == "POLR2B"
  assert "NOGENE_HUMAN" not in mappings
