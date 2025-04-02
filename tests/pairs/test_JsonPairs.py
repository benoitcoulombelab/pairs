import json
import os
import re
import shutil
from pathlib import Path
from unittest.mock import MagicMock, ANY

import pytest
from Bio import SeqIO

from pairs import JsonPairs


@pytest.fixture
def mock_testclass():
  _json_pairs = JsonPairs.json_pairs
  _parse_fasta = JsonPairs.parse_fasta
  yield
  JsonPairs.json_pairs = _json_pairs
  JsonPairs.parse_fasta = _parse_fasta


def test_main(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  JsonPairs.json_pairs = MagicMock()
  JsonPairs.main(["-b", baits, "-t", targets])
  JsonPairs.json_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                               seeds=None,
                                               unique=False,
                                               skip_identity=False,
                                               output="")
  assert JsonPairs.json_pairs.call_args.kwargs["baits"].name == baits
  assert JsonPairs.json_pairs.call_args.kwargs["targets"].name == targets


def test_main_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  seed1 = 12
  seed2 = 30
  output = "output_dir"
  os.mkdir(output)
  JsonPairs.json_pairs = MagicMock()
  JsonPairs.main(
      ["-b", baits, "-t", targets, "-s", str(seed1), str(seed2), "-u", "-i",
       "-o", output])
  JsonPairs.json_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                               seeds=[seed1, seed2],
                                               unique=True,
                                               skip_identity=True,
                                               output=output)
  assert JsonPairs.json_pairs.call_args.kwargs["baits"].name == baits
  assert JsonPairs.json_pairs.call_args.kwargs["targets"].name == targets


def test_main_long_parameters(testdir, mock_testclass):
  baits = "baits.fasta"
  open(baits, 'w').close()
  targets = "targets.fasta"
  open(targets, 'w').close()
  seed1 = 12
  seed2 = 30
  output = "output_dir"
  os.mkdir(output)
  JsonPairs.json_pairs = MagicMock()
  JsonPairs.main(
      ["--baits", baits, "--targets", targets, "--seed", str(seed1), str(seed2),
       "--unique", "--identity", "--output", output])
  JsonPairs.json_pairs.assert_called_once_with(baits=ANY, targets=ANY,
                                               seeds=[seed1, seed2],
                                               unique=True,
                                               skip_identity=True,
                                               output=output)
  assert JsonPairs.json_pairs.call_args.kwargs["baits"].name == baits
  assert JsonPairs.json_pairs.call_args.kwargs["targets"].name == targets


def test_json_pairs(testdir, mock_testclass):
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
  JsonPairs.json_pairs(baits=baits, targets=targets)
  assert os.path.isfile("RPAB1_HUMAN__RPB4_HUMAN.json")
  with open("RPAB1_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPAB1_HUMAN__RPB7_HUMAN.json")
  with open("RPAB1_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq
  assert os.path.isfile("RPB9_HUMAN__RPB4_HUMAN.json")
  with open("RPB9_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPB9_HUMAN__RPB7_HUMAN.json")
  with open("RPB9_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq


def test_json_pairs_seed(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  seed = 12
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
  JsonPairs.json_pairs(baits=baits, targets=targets, seeds=[seed])
  assert os.path.isfile("RPAB1_HUMAN__RPB4_HUMAN.json")
  with open("RPAB1_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert json_data["modelSeeds"][0] == seed
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPAB1_HUMAN__RPB7_HUMAN.json")
  with open("RPAB1_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert json_data["modelSeeds"][0] == seed
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq
  assert os.path.isfile("RPB9_HUMAN__RPB4_HUMAN.json")
  with open("RPB9_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert json_data["modelSeeds"][0] == seed
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPB9_HUMAN__RPB7_HUMAN.json")
  with open("RPB9_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert json_data["modelSeeds"][0] == seed
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq


def test_json_pairs_seeds(testdir, mock_testclass):
  baits = "baits.fasta"
  targets = "targets.fasta"
  seed1 = 12
  seed2 = 30
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
  JsonPairs.json_pairs(baits=baits, targets=targets, seeds=[seed1, seed2])
  assert os.path.isfile("RPAB1_HUMAN__RPB4_HUMAN.json")
  with open("RPAB1_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 2
    assert json_data["modelSeeds"][0] == seed1
    assert json_data["modelSeeds"][1] == seed2
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPAB1_HUMAN__RPB7_HUMAN.json")
  with open("RPAB1_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 2
    assert json_data["modelSeeds"][0] == seed1
    assert json_data["modelSeeds"][1] == seed2
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq
  assert os.path.isfile("RPB9_HUMAN__RPB4_HUMAN.json")
  with open("RPB9_HUMAN__RPB4_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 2
    assert json_data["modelSeeds"][0] == seed1
    assert json_data["modelSeeds"][1] == seed2
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile("RPB9_HUMAN__RPB7_HUMAN.json")
  with open("RPB9_HUMAN__RPB7_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 2
    assert json_data["modelSeeds"][0] == seed1
    assert json_data["modelSeeds"][1] == seed2
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq


def test_json_pairs_unique(testdir, mock_testclass):
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
  JsonPairs.json_pairs(baits=baits, targets=targets, unique=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.json")
  with open("RPAB1_HUMAN__RPB9_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB9_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2i.seq
  assert os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.json")
  with open("RPAB1_HUMAN__RPAB1_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPAB1_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPABA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2e.seq
  assert not os.path.isfile("RPB9_HUMAN__RPAB1_HUMAN.json")
  assert os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.json")
  with open("RPB9_HUMAN__RPB9_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB9_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2i.seq


def test_json_pairs_skip_identity(testdir, mock_testclass):
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
  JsonPairs.json_pairs(baits=baits, targets=targets, skip_identity=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.json")
  with open("RPAB1_HUMAN__RPB9_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB9_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2i.seq
  assert not os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.json")
  assert os.path.isfile("RPB9_HUMAN__RPAB1_HUMAN.json")
  with open("RPB9_HUMAN__RPAB1_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPAB1_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2e.seq
  assert not os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.json")


def test_json_pairs_unique_skip_identity(testdir, mock_testclass):
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
  JsonPairs.json_pairs(baits=baits, targets=targets, unique=True,
                       skip_identity=True)
  assert os.path.isfile("RPAB1_HUMAN__RPB9_HUMAN.json")
  with open("RPAB1_HUMAN__RPB9_HUMAN.json", 'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB9_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2i.seq
  assert not os.path.isfile("RPAB1_HUMAN__RPAB1_HUMAN.json")
  assert not os.path.isfile("RPB9_HUMAN__RPAB1_HUMAN.json")
  assert not os.path.isfile("RPB9_HUMAN__RPB9_HUMAN.json")


def test_json_pairs_output(testdir, mock_testclass):
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
  JsonPairs.json_pairs(baits=baits, targets=targets, output=output)
  assert os.path.isfile(os.path.join(output, "RPAB1_HUMAN__RPB4_HUMAN.json"))
  with open(os.path.join(output, "RPAB1_HUMAN__RPB4_HUMAN.json"),
            'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile(os.path.join(output, "RPAB1_HUMAN__RPB7_HUMAN.json"))
  with open(os.path.join(output, "RPAB1_HUMAN__RPB7_HUMAN.json"),
            'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPAB1_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPAB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2e.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq
  assert os.path.isfile(os.path.join(output, "RPB9_HUMAN__RPB4_HUMAN.json"))
  with open(os.path.join(output, "RPB9_HUMAN__RPB4_HUMAN.json"),
            'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB4_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2d.seq
  assert os.path.isfile(os.path.join(output, "RPB9_HUMAN__RPB7_HUMAN.json"))
  with open(os.path.join(output, "RPB9_HUMAN__RPB7_HUMAN.json"),
            'r') as json_file:
    json_data = json.load(json_file)
    assert "name" in json_data
    assert json_data["name"] == "RPB9_HUMAN__RPB7_HUMAN"
    assert "modelSeeds" in json_data
    assert len(json_data["modelSeeds"]) == 1
    assert re.match(r"^\d+$", str(json_data["modelSeeds"][0]))
    assert "dialect" in json_data
    assert json_data["dialect"] == "alphafold3"
    assert "version" in json_data
    assert json_data["version"] == 1
    assert "sequences" in json_data
    assert len(json_data["sequences"]) == 2
    assert "protein" in json_data["sequences"][0]
    assert "id" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["id"] == "RPB"
    assert "sequence" in json_data["sequences"][0]["protein"]
    assert json_data["sequences"][0]["protein"]["sequence"] == polr2i.seq
    assert "protein" in json_data["sequences"][1]
    assert "id" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["id"] == "RPBA"
    assert "sequence" in json_data["sequences"][1]["protein"]
    assert json_data["sequences"][1]["protein"]["sequence"] == polr2g.seq
