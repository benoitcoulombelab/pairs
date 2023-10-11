import os.path
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import pytest
from Bio.Seq import Seq

from afpairs import SplitFasta, FastaId


@pytest.fixture
def mock_testclass():
    _split_fasta = SplitFasta.split_fasta
    _parse_fasta = SplitFasta.parse_fasta
    _dir_path = SplitFasta.dir_path
    _fasta_id = FastaId.fasta_id
    yield
    SplitFasta.split_fasta = _split_fasta
    SplitFasta.parse_fasta = _parse_fasta
    SplitFasta.dir_path = _dir_path
    FastaId.fasta_id = _fasta_id


def test_main(testdir, mock_testclass):
    SplitFasta.split_fasta = MagicMock()
    stdin_file = "stdin.txt"
    open(stdin_file, 'w').close()
    with open(stdin_file, 'r') as stdin_in, patch('sys.stdin', stdin_in):
        SplitFasta.main()
    SplitFasta.split_fasta.assert_called_once_with(
        fasta=ANY, output_dir="")
    fasta = SplitFasta.split_fasta.call_args.kwargs["fasta"]
    assert isinstance(fasta, TextIOWrapper)
    assert fasta.mode == "r"


def test_main_parameters(testdir, mock_testclass):
    fasta = "abc.fasta"
    open(fasta, 'w').close()
    output_dir = "output"
    os.mkdir(output_dir)
    SplitFasta.split_fasta = MagicMock()
    SplitFasta.main(["-o", output_dir, fasta])
    SplitFasta.split_fasta.assert_called_once_with(
        fasta=ANY, output_dir=output_dir)
    fasta_in = SplitFasta.split_fasta.call_args.kwargs["fasta"]
    assert fasta_in.name == fasta
    assert fasta_in.mode == "r"


def test_main_long_parameters(testdir, mock_testclass):
    fasta = "abc.fasta"
    open(fasta, 'w').close()
    output_dir = "output"
    os.mkdir(output_dir)
    SplitFasta.split_fasta = MagicMock()
    SplitFasta.main(["--output", output_dir, fasta])
    SplitFasta.split_fasta.assert_called_once_with(
        fasta=ANY, output_dir=output_dir)
    fasta_in = SplitFasta.split_fasta.call_args.kwargs["fasta"]
    assert fasta_in.name == fasta
    assert fasta_in.mode == "r"


def test_split_fasta(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    output_dir = "output"
    os.mkdir(output_dir)
    SplitFasta.split_fasta(fasta=fasta, output_dir=output_dir)
    assert os.path.isfile(f"{output_dir}/RPAB1_HUMAN.fasta")
    with open(f"{output_dir}/RPAB1_HUMAN.fasta", 'r') as fasta_in:
        assert fasta_in.readline() == ">sp|P19388|RPAB1_HUMAN " \
                                      "DNA-directed RNA polymerases I, II, and III subunit RPABC1 " \
                                      "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4\n"
        assert fasta_in.readline() == "MDDEEETYRLWKIRKTIMQLCHDRGYLVTQDELDQTLEEFKAQSGDKPSEGRPRRTDLTV\n"
        assert fasta_in.readline() == "LVAHNDDPTDQMFVFFPEEPKVGIKTIKVYCQRMQEENITRALIVVQQGMTPSAKQSLVD\n"
        assert fasta_in.readline() == "MAPKYILEQFLQQELLINITEHELVPEHVVMTKEEVTELLARYKLRENQLPRIQAGDPVA\n"
        assert fasta_in.readline() == "RYFGIKRGQVVKIIRPSETAGRYITYRLVQ\n"
    assert os.path.isfile(f"{output_dir}/RPB9_HUMAN.fasta")
    with open(f"{output_dir}/RPB9_HUMAN.fasta", 'r') as fasta_in:
        assert fasta_in.readline() == ">sp|P36954|RPB9_HUMAN " \
                                      "DNA-directed RNA polymerase II subunit RPB9 " \
                                      "OS=Homo sapiens OX=9606 GN=POLR2I PE=1 SV=1\n"
        assert fasta_in.readline() == "MEPDGTYEPGFVGIRFCQECNNMLYPKEDKENRILLYACRNCDYQQEADNSCIYVNKITH\n"
        assert fasta_in.readline() == "EVDELTQIIADVSQDPTLPRTEDHPCQKCGHKEAVFFQSHSARAEDAMRLYYVCTAPHCG\n"
        assert fasta_in.readline() == "HRWTE\n"


def test_parse_fasta(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    FastaId.fasta_id = MagicMock(side_effect=["RPB1_HUMAN", "RPB2_HUMAN"])
    sequences = SplitFasta.parse_fasta(fasta)
    FastaId.fasta_id.assert_any_call(
        "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
        "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4")
    FastaId.fasta_id.assert_any_call(
        "sp|P36954|RPB9_HUMAN DNA-directed RNA polymerase II subunit RPB9 "
        "OS=Homo sapiens OX=9606 GN=POLR2I PE=1 SV=1")
    assert len(sequences) == 2
    assert "RPB1_HUMAN" in sequences
    assert sequences["RPB1_HUMAN"].name == "sp|P19388|RPAB1_HUMAN"
    assert sequences["RPB1_HUMAN"].seq == Seq("MDDEEETYRLWKIRKTIMQLCHDRGYLVTQDELDQTLEEFKAQSGDKPSEGRPRRTDLTV"
                                              "LVAHNDDPTDQMFVFFPEEPKVGIKTIKVYCQRMQEENITRALIVVQQGMTPSAKQSLVD"
                                              "MAPKYILEQFLQQELLINITEHELVPEHVVMTKEEVTELLARYKLRENQLPRIQAGDPVA"
                                              "RYFGIKRGQVVKIIRPSETAGRYITYRLVQ")
    assert "RPB2_HUMAN" in sequences
    assert sequences["RPB2_HUMAN"].name == "sp|P36954|RPB9_HUMAN"
    assert sequences["RPB2_HUMAN"].seq == Seq("MEPDGTYEPGFVGIRFCQECNNMLYPKEDKENRILLYACRNCDYQQEADNSCIYVNKITH"
                                              "EVDELTQIIADVSQDPTLPRTEDHPCQKCGHKEAVFFQSHSARAEDAMRLYYVCTAPHCG"
                                              "HRWTE")
