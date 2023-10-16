from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pairs import FastaId


@pytest.fixture
def mock_testclass():
    _fasta_id = FastaId.fasta_id
    yield
    FastaId.fasta_id = _fasta_id


def test_main(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    FastaId.fasta_id = MagicMock(return_value="P19388")
    FastaId.main([str(fasta)])
    FastaId.fasta_id.assert_called_once_with(
        "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
        "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4",
        gene=False)


def test_main_2(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P36954.fasta")
    FastaId.fasta_id = MagicMock(return_value="P36954")
    FastaId.main([str(fasta)])
    FastaId.fasta_id.assert_called_once_with(
        "sp|P36954|RPB9_HUMAN DNA-directed RNA polymerase II subunit RPB9 "
        "OS=Homo sapiens OX=9606 GN=POLR2I PE=1 SV=1",
        gene=False)


def test_main_no_name(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("no_name.fasta")
    FastaId.fasta_id = MagicMock(return_value="P19388")
    try:
        FastaId.main([str(fasta)])
        assert False, "Should raise AssertionError"
    except AssertionError:
        pass
    finally:
        FastaId.fasta_id.assert_not_called()


def test_main_gene(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    FastaId.fasta_id = MagicMock(return_value="P19388")
    FastaId.main([str(fasta), "-g"])
    FastaId.fasta_id.assert_called_once_with(
        "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
        "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4",
        gene=True)


def test_main_longgene(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    FastaId.fasta_id = MagicMock(return_value="P19388")
    FastaId.main([str(fasta), "--gene"])
    FastaId.fasta_id.assert_called_once_with(
        "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
        "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4",
        gene=True)


def test_fasta_id_uniprot_description(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
                                "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4")
    assert fasta_id == "RPAB1_HUMAN"


def test_fasta_id_uniprot_description_no_greater(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id("sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
                                "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4")
    assert fasta_id == "RPAB1_HUMAN"


def test_fasta_id_uniprot_description_gene(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
                                "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4", gene=True)
    assert fasta_id == "POLR2E"


def test_fasta_id_uniprot(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">sp|P19388|RPAB1_HUMAN")
    assert fasta_id == "RPAB1_HUMAN"


def test_fasta_id_uniprot_no_greater(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id("sp|P19388|RPAB1_HUMAN")
    assert fasta_id == "RPAB1_HUMAN"


def test_fasta_id_uniprot_gene(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id("sp|P19388|RPAB1_HUMAN", gene=True)
    assert fasta_id == "RPAB1_HUMAN"


def test_fasta_id_uniprot_hyphen(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id("sp|P19388|RPAB1-HUMAN")
    assert fasta_id == "RPAB1-HUMAN"


def test_fasta_id_refseq(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">NP_000928.1")
    assert fasta_id == "NP_000928"


def test_fasta_id_refseq_gene(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">NP_000928.1", gene=True)
    assert fasta_id == "NP_000928"


def test_fasta_id_refseq_gi(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id(">gi|4505939|ref|NP_000928.1|")
    assert fasta_id == "NP_000928"


def test_fasta_id_none(testdir, mock_testclass):
    fasta_id = FastaId.fasta_id("")
    assert fasta_id is None
