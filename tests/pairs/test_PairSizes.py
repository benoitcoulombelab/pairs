import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock, ANY

import pytest

from pairs import PairSizes


@pytest.fixture
def mock_testclass():
  _pair_sizes = PairSizes.pair_sizes
  _parse_fasta = PairSizes.parse_fasta
  yield
  PairSizes.pair_sizes = _pair_sizes
  PairSizes.parse_fasta = _parse_fasta


def test_main(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  PairSizes.pair_sizes = MagicMock()
  PairSizes.main(["-b", baits, "-t", targets])
  PairSizes.pair_sizes.assert_called_once_with(baits=ANY, targets=ANY,
                                               unique=False,
                                               skip_identity=False,
                                               output=ANY)
  assert PairSizes.pair_sizes.call_args.kwargs["baits"].name == baits
  assert PairSizes.pair_sizes.call_args.kwargs["targets"].name == targets
  assert PairSizes.pair_sizes.call_args.kwargs[
           "output"].name == "pair_sizes.txt"


def test_main_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  sizes = "sizes.txt"
  PairSizes.pair_sizes = MagicMock()
  PairSizes.main(
      ["-b", baits, "-t", targets, "-u", "-i",
       "-o", sizes])
  PairSizes.pair_sizes.assert_called_once_with(baits=ANY, targets=ANY,
                                               unique=True,
                                               skip_identity=True,
                                               output=ANY)
  assert PairSizes.pair_sizes.call_args.kwargs["baits"].name == baits
  assert PairSizes.pair_sizes.call_args.kwargs["targets"].name == targets
  assert PairSizes.pair_sizes.call_args.kwargs["output"].name == sizes


def test_main_long_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  sizes = "sizes.txt"
  PairSizes.pair_sizes = MagicMock()
  PairSizes.main(
      ["--baits", baits, "--targets", targets,
       "--unique", "--identity", "--output", sizes])
  PairSizes.pair_sizes.assert_called_once_with(baits=ANY, targets=ANY,
                                               unique=True,
                                               skip_identity=True,
                                               output=ANY)
  assert PairSizes.pair_sizes.call_args.kwargs["baits"].name == baits
  assert PairSizes.pair_sizes.call_args.kwargs["targets"].name == targets
  assert PairSizes.pair_sizes.call_args.kwargs["output"].name == sizes


def test_pair_sizes(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2d_file = Path(__file__).parent.joinpath("O15514.fasta")
  polr2g_file = Path(__file__).parent.joinpath("P62487.fasta")
  with open(baits, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  with open(targets, "wb") as output:
    with open(polr2d_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2g_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  sizes = "pair_sizes.txt"
  with open(sizes, 'w') as sizes_file:
    PairSizes.pair_sizes(baits=baits, targets=targets, output=sizes_file)
  assert os.path.isfile(sizes)
  with open(sizes, 'r') as sizes_file:
    assert sizes_file.readline() == "RPAB1_HUMAN__RPB4_HUMAN\t352\n"
    assert sizes_file.readline() == "RPAB1_HUMAN__RPB7_HUMAN\t382\n"
    assert sizes_file.readline() == "RPB9_HUMAN__RPB4_HUMAN\t267\n"
    assert sizes_file.readline() == "RPB9_HUMAN__RPB7_HUMAN\t297\n"
    assert sizes_file.readline() == ""


def test_pair_sizes_unique(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  with open(baits, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  with open(targets, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  sizes = "pair_sizes.txt"
  with open(sizes, 'w') as sizes_file:
    PairSizes.pair_sizes(baits=baits, targets=targets, output=sizes_file,
                         unique=True)
  assert os.path.isfile(sizes)
  with open(sizes, 'r') as sizes_file:
    assert sizes_file.readline() == "RPAB1_HUMAN__RPAB1_HUMAN\t420\n"
    assert sizes_file.readline() == "RPAB1_HUMAN__RPB9_HUMAN\t335\n"
    assert sizes_file.readline() == "RPB9_HUMAN__RPB9_HUMAN\t250\n"
    assert sizes_file.readline() == ""


def test_pair_sizes_skip_identity(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  with open(baits, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  with open(targets, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  sizes = "pair_sizes.txt"
  with open(sizes, 'w') as sizes_file:
    PairSizes.pair_sizes(baits=baits, targets=targets, output=sizes_file,
                         skip_identity=True)
  assert os.path.isfile(sizes)
  with open(sizes, 'r') as sizes_file:
    assert sizes_file.readline() == "RPAB1_HUMAN__RPB9_HUMAN\t335\n"
    assert sizes_file.readline() == "RPB9_HUMAN__RPAB1_HUMAN\t335\n"
    assert sizes_file.readline() == ""


def test_pair_sizes_unique_skip_identity(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  with open(baits, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  with open(targets, "wb") as output:
    with open(polr2e_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
    with open(polr2i_file, "rb") as infile:
      shutil.copyfileobj(infile, output)
  sizes = "pair_sizes.txt"
  with open(sizes, 'w') as sizes_file:
    PairSizes.pair_sizes(baits=baits, targets=targets, output=sizes_file,
                         unique=True, skip_identity=True)
  assert os.path.isfile(sizes)
  with open(sizes, 'r') as sizes_file:
    assert sizes_file.readline() == "RPAB1_HUMAN__RPB9_HUMAN\t335\n"
    assert sizes_file.readline() == ""
