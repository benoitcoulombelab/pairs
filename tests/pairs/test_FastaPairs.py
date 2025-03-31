import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock, ANY

import pytest
from Bio import SeqIO

from pairs import FastaPairs


@pytest.fixture
def mock_testclass():
  _fasta_pairs = FastaPairs.fasta_pairs
  _parse_fasta = FastaPairs.parse_fasta
  yield
  FastaPairs.fasta_pairs = _fasta_pairs
  FastaPairs.parse_fasta = _parse_fasta


def test_main(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  FastaPairs.fasta_pairs = MagicMock()
  FastaPairs.main(["-b", baits, "-t", targets])
  FastaPairs.fasta_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                                 unique=False,
                                                 skip_identity=False,
                                                 output="")
  assert FastaPairs.fasta_pairs.call_args.kwargs["baits"].name == baits
  assert FastaPairs.fasta_pairs.call_args.kwargs["targets"].name == targets


def test_main_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  output = "output_dir"
  os.mkdir(output)
  FastaPairs.fasta_pairs = MagicMock()
  FastaPairs.main(["-b", baits, "-t", targets, "-u", "-i", "-o", output])
  FastaPairs.fasta_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                                 unique=True,
                                                 skip_identity=True,
                                                 output=output)
  assert FastaPairs.fasta_pairs.call_args.kwargs["baits"].name == baits
  assert FastaPairs.fasta_pairs.call_args.kwargs["targets"].name == targets


def test_main_long_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  output = "output_dir"
  os.mkdir(output)
  FastaPairs.fasta_pairs = MagicMock()
  FastaPairs.main(
      ["--baits", baits, "--targets", targets, "--unique", "--identity",
       "--output", output])
  FastaPairs.fasta_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                                 unique=True,
                                                 skip_identity=True,
                                                 output=output)
  assert FastaPairs.fasta_pairs.call_args.kwargs["baits"].name == baits
  assert FastaPairs.fasta_pairs.call_args.kwargs["targets"].name == targets


def test_fasta_pairs(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2e = next(SeqIO.parse(polr2e_file, "fasta"))
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2i = next(SeqIO.parse(polr2i_file, "fasta"))
  polr2d_file = Path(__file__).parent.joinpath("O15514.fasta")
  polr2d = next(SeqIO.parse(polr2d_file, "fasta"))
  polr2g_file = Path(__file__).parent.joinpath("P62487.fasta")
  polr2g = next(SeqIO.parse(polr2g_file, "fasta"))
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
  FastaPairs.fasta_pairs(baits=baits, targets=targets)
  assert os.path.isfile("RPAB1_HUMAN__RPB4_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPB4_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2d.description
  assert sequences[1].seq == polr2d.seq
  assert os.path.isfile("RPAB1_HUMAN__RPB7_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPB7_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2g.description
  assert sequences[1].seq == polr2g.seq
  assert os.path.isfile("RPB9_HUMAN__RPB4_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPB9_HUMAN__RPB4_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2d.description
  assert sequences[1].seq == polr2d.seq
  assert os.path.isfile("RPB9_HUMAN__RPB7_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPB9_HUMAN__RPB7_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2g.description
  assert sequences[1].seq == polr2g.seq


def test_fasta_pairs_unique(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2e = next(SeqIO.parse(polr2e_file, "fasta"))
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2i = next(SeqIO.parse(polr2i_file, "fasta"))
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
  FastaPairs.fasta_pairs(baits=baits, targets=targets, unique=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPB9_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2i.description
  assert sequences[1].seq == polr2i.seq
  assert os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPAB1_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2e.description
  assert sequences[1].seq == polr2e.seq
  assert not os.path.isfile("RPB9_HUMAN__RPAB1_HUMAN.fasta")
  assert os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPB9_HUMAN__RPB9_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2i.description
  assert sequences[1].seq == polr2i.seq


def test_fasta_pairs_skip_identity(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2e = next(SeqIO.parse(polr2e_file, "fasta"))
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2i = next(SeqIO.parse(polr2i_file, "fasta"))
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
  FastaPairs.fasta_pairs(baits=baits, targets=targets, skip_identity=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPB9_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2i.description
  assert sequences[1].seq == polr2i.seq
  assert not os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPB9_HUMAN__RPAB1_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2e.description
  assert sequences[1].seq == polr2e.seq
  assert not os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.fasta")


def test_fasta_pairs_unique_skip_identity(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2e = next(SeqIO.parse(polr2e_file, "fasta"))
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2i = next(SeqIO.parse(polr2i_file, "fasta"))
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
  FastaPairs.fasta_pairs(baits=baits, targets=targets, unique=True,
                         skip_identity=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.fasta")
  sequences = list(SeqIO.parse("RPAB1_HUMAN__RPB9_HUMAN.fasta", "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2i.description
  assert sequences[1].seq == polr2i.seq
  assert not os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.fasta")
  assert not os.path.isfile("RPB9_HUMAN__RPAB1_HUMAN.fasta")
  assert not os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.fasta")


def test_fasta_pairs_output(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  polr2e_file = Path(__file__).parent.joinpath("P19388.fasta")
  polr2e = next(SeqIO.parse(polr2e_file, "fasta"))
  polr2i_file = Path(__file__).parent.joinpath("P36954.fasta")
  polr2i = next(SeqIO.parse(polr2i_file, "fasta"))
  polr2d_file = Path(__file__).parent.joinpath("O15514.fasta")
  polr2d = next(SeqIO.parse(polr2d_file, "fasta"))
  polr2g_file = Path(__file__).parent.joinpath("P62487.fasta")
  polr2g = next(SeqIO.parse(polr2g_file, "fasta"))
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
  output = "output_dir"
  os.mkdir(output)
  FastaPairs.fasta_pairs(baits=baits, targets=targets, output=output)
  assert os.path.isfile(os.path.join(output, "RPAB1_HUMAN__RPB4_HUMAN.fasta"))
  sequences = list(
      SeqIO.parse(os.path.join(output, "RPAB1_HUMAN__RPB4_HUMAN.fasta"),
                  "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2d.description
  assert sequences[1].seq == polr2d.seq
  assert os.path.isfile(os.path.join(output, "RPAB1_HUMAN__RPB7_HUMAN.fasta"))
  sequences = list(
      SeqIO.parse(os.path.join(output, "RPAB1_HUMAN__RPB7_HUMAN.fasta"),
                  "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2e.description
  assert sequences[0].seq == polr2e.seq
  assert sequences[1].description == polr2g.description
  assert sequences[1].seq == polr2g.seq
  assert os.path.isfile(os.path.join(output, "RPB9_HUMAN__RPB4_HUMAN.fasta"))
  sequences = list(
      SeqIO.parse(os.path.join(output, "RPB9_HUMAN__RPB4_HUMAN.fasta"),
                  "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2d.description
  assert sequences[1].seq == polr2d.seq
  assert os.path.isfile(os.path.join(output, "RPB9_HUMAN__RPB7_HUMAN.fasta"))
  sequences = list(
      SeqIO.parse(os.path.join(output, "RPB9_HUMAN__RPB7_HUMAN.fasta"),
                  "fasta"))
  assert len(sequences) == 2
  assert sequences[0].description == polr2i.description
  assert sequences[0].seq == polr2i.seq
  assert sequences[1].description == polr2g.description
  assert sequences[1].seq == polr2g.seq
