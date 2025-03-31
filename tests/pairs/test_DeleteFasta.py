import os.path
from pathlib import Path
from shutil import copyfile
from unittest.mock import MagicMock

import pytest

from pairs import DeleteFasta


@pytest.fixture
def mock_testclass():
  _delete_fasta = DeleteFasta.delete_fasta
  _parse_fasta = DeleteFasta.parse_fasta
  yield
  DeleteFasta.delete_fasta = _delete_fasta
  DeleteFasta.parse_fasta = _parse_fasta


def test_main(testdir, mock_testclass):
  DeleteFasta.delete_fasta = MagicMock()
  DeleteFasta.main(["P19388__P36954.fasta", "P19388.fasta"])
  DeleteFasta.delete_fasta.assert_called_once_with(
      inputs=["P19388__P36954.fasta", "P19388.fasta"],
      invalid_sequence=False,
      length=None, backup=None, verbose=False)


def test_main_parameters(testdir, mock_testclass):
  os.mkdir("backup")
  DeleteFasta.delete_fasta = MagicMock()
  DeleteFasta.main(
      ["-s", "-l", "200", "-b", "backup", "-v", "P19388__P36954.fasta",
       "P19388.fasta"])
  DeleteFasta.delete_fasta.assert_called_once_with(
      inputs=["P19388__P36954.fasta", "P19388.fasta"],
      invalid_sequence=True,
      length=200, backup="backup", verbose=True)


def test_main_long_parameters(testdir, mock_testclass):
  os.mkdir("backup")
  DeleteFasta.delete_fasta = MagicMock()
  DeleteFasta.main(
      ["--sequence", "--length", "200", "--backup", "backup", "--verbose",
       "P19388__P36954.fasta", "P19388.fasta"])
  DeleteFasta.delete_fasta.assert_called_once_with(
      inputs=["P19388__P36954.fasta", "P19388.fasta"],
      invalid_sequence=True,
      length=200, backup="backup", verbose=True)


def test_delete_fasta_no_action(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  copyfile(fasta1, "P19388__P36954.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  DeleteFasta.delete_fasta(["P19388__P36954.fasta", "P19388.fasta"])
  assert os.path.isfile("P19388__P36954.fasta")
  assert os.path.isfile("P19388.fasta")


def test_delete_fasta_delete_invalid_sequence(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P01786.fasta")
  copyfile(fasta1, "P01786.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  fasta3 = Path(__file__).parent.joinpath("P01788.fasta")
  copyfile(fasta3, "P01788.fasta")
  DeleteFasta.delete_fasta(["P01786.fasta", "P19388.fasta", "P01788.fasta"],
                           invalid_sequence=True)
  assert not os.path.isfile("P01786.fasta")
  assert os.path.isfile("P19388.fasta")
  assert not os.path.isfile("P01788.fasta")


def test_delete_fasta_delete_length_all(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  copyfile(fasta1, "P19388__P36954.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  DeleteFasta.delete_fasta(["P19388__P36954.fasta", "P19388.fasta"], length=10)
  assert not os.path.isfile("P19388__P36954.fasta")
  assert not os.path.isfile("P19388.fasta")


def test_delete_fasta_delete_length_one(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  copyfile(fasta1, "P19388__P36954.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  DeleteFasta.delete_fasta(["P19388__P36954.fasta", "P19388.fasta"], length=250)
  assert not os.path.isfile("P19388__P36954.fasta")
  assert os.path.isfile("P19388.fasta")


def test_delete_fasta_delete_length_zero(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  copyfile(fasta1, "P19388__P36954.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  DeleteFasta.delete_fasta(["P19388__P36954.fasta", "P19388.fasta"], length=0)
  assert not os.path.isfile("P19388__P36954.fasta")
  assert not os.path.isfile("P19388.fasta")


def test_delete_fasta_backup(testdir, mock_testclass):
  fasta1 = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  copyfile(fasta1, "P19388__P36954.fasta")
  fasta2 = Path(__file__).parent.joinpath("P19388.fasta")
  copyfile(fasta2, "P19388.fasta")
  os.mkdir("backup")
  DeleteFasta.delete_fasta(["P19388__P36954.fasta", "P19388.fasta"], length=250,
                           backup="backup")
  assert not os.path.isfile("P19388__P36954.fasta")
  assert os.path.isfile("P19388.fasta")
  assert os.path.isfile("backup/P19388__P36954.fasta")
  assert not os.path.isfile("backup/P19388.fasta")


def test_parse_fasta(testdir, mock_testclass):
  fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
  sequences = DeleteFasta.parse_fasta(fasta)
  assert sequences[0].name == "sp|P19388|RPAB1_HUMAN"
  assert sequences[0].description == (
    "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
    "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4")
  assert sequences[0].seq == (
    "MDDEEETYRLWKIRKTIMQLCHDRGYLVTQDELDQTLEEFKAQSGDKPSEGRPRRTDLTV"
    "LVAHNDDPTDQMFVFFPEEPKVGIKTIKVYCQRMQEENITRALIVVQQGMTPSAKQSLVD"
    "MAPKYILEQFLQQELLINITEHELVPEHVVMTKEEVTELLARYKLRENQLPRIQAGDPVA"
    "RYFGIKRGQVVKIIRPSETAGRYITYRLVQ")
  assert sequences[1].name == "sp|P36954|RPB9_HUMAN"
  assert sequences[1].description == (
    "sp|P36954|RPB9_HUMAN DNA-directed RNA polymerase II subunit RPB9 "
    "OS=Homo sapiens OX=9606 GN=POLR2I PE=1 SV=1")
  assert sequences[1].seq == (
    "MEPDGTYEPGFVGIRFCQECNNMLYPKEDKENRILLYACRNCDYQQEADNSCIYVNKITH"
    "EVDELTQIIADVSQDPTLPRTEDHPCQKCGHKEAVFFQSHSARAEDAMRLYYVCTAPHCG"
    "HRWTE")
