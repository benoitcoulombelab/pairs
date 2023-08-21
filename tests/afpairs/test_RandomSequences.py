import os.path
import sys
from pathlib import Path
from unittest.mock import MagicMock, ANY

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from afpairs import RandomSequences


@pytest.fixture
def mock_testclass():
    _random_sequences = RandomSequences.random_sequences
    _generate_sequence = RandomSequences.generate_sequence
    _aa_count = RandomSequences.aa_count
    _parse_fasta = RandomSequences.parse_fasta
    yield
    RandomSequences.random_sequences = _random_sequences
    RandomSequences.generate_sequence = _generate_sequence
    RandomSequences.aa_count = _aa_count
    RandomSequences.parse_fasta = _parse_fasta


def test_main(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    RandomSequences.random_sequences = MagicMock()
    RandomSequences.main([str(fasta)])
    RandomSequences.random_sequences.assert_called_once_with(
        fasta=ANY, output=sys.stdout, invalid_sequence=False)
    assert RandomSequences.random_sequences.call_args.kwargs["fasta"].name == str(fasta)


def test_main_parameters(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    RandomSequences.random_sequences = MagicMock()
    RandomSequences.main(["-s", str(fasta), "output.fasta"])
    RandomSequences.random_sequences.assert_called_once_with(
        fasta=ANY, output=ANY, invalid_sequence=True)
    assert RandomSequences.random_sequences.call_args.kwargs["fasta"].name == str(fasta)
    assert RandomSequences.random_sequences.call_args.kwargs["output"].name == "output.fasta"


def test_main_long_parameters(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    RandomSequences.random_sequences = MagicMock()
    RandomSequences.main(["--sequence", str(fasta), "output.fasta"])
    RandomSequences.random_sequences.assert_called_once_with(
        fasta=ANY, output=ANY, invalid_sequence=True)
    assert RandomSequences.random_sequences.call_args.kwargs["fasta"].name == str(fasta)
    assert RandomSequences.random_sequences.call_args.kwargs["output"].name == "output.fasta"


def test_random_sequences(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    output = "output.fasta"
    sequences = [SeqRecord(Seq("MFYHISLEHEILLHPRYF"), id="P62487"),
                 SeqRecord(Seq("MEPDGTYEPGFVGIRFCQECNNML"), id="P36954")]
    aa_count = {'A': 17, 'C': 11, 'D': 21, 'E': 33}
    RandomSequences.parse_fasta = MagicMock(return_value=sequences)
    RandomSequences.aa_count = MagicMock(return_value=aa_count)
    RandomSequences.generate_sequence = MagicMock(
        side_effect=[Seq("EVKLVESGGGLVQPGGSL"), Seq("ALEWLGFIRNKAOGYTT")])
    with open(output, "w") as output_file:
        RandomSequences.random_sequences(fasta, output_file)
    RandomSequences.parse_fasta.assert_called_once_with(fasta)
    RandomSequences.aa_count.assert_called_once_with(sequences)
    RandomSequences.generate_sequence.assert_called()
    assert RandomSequences.generate_sequence.call_count == 2
    assert RandomSequences.generate_sequence.call_args_list[0].args[0] == len(sequences[0].seq)
    aa_probabilities = RandomSequences.generate_sequence.call_args_list[0].args[1]
    assert len(aa_probabilities) == 4
    assert aa_probabilities["A"] - 0.207 < 0.001
    assert aa_probabilities["C"] - 0.134 < 0.001
    assert aa_probabilities["D"] - 0.256 < 0.001
    assert aa_probabilities["E"] - 0.402 < 0.001
    assert RandomSequences.generate_sequence.call_args_list[1].args[0] == len(sequences[1].seq)
    assert RandomSequences.generate_sequence.call_args_list[1].args[1] == \
           RandomSequences.generate_sequence.call_args_list[0].args[1]
    assert os.path.isfile(output)
    out_sequences = []
    for record in SeqIO.parse(output, "fasta"):
        out_sequences.append(record)
    with open(output, "r") as output_file:
        print(output_file.readlines())
    assert len(out_sequences) == 2
    assert out_sequences[0].name == "DECOY_1"
    assert out_sequences[0].description == "DECOY_1 MASS=1725"
    assert out_sequences[0].seq == "EVKLVESGGGLVQPGGSL"
    assert out_sequences[1].name == "DECOY_2"
    assert out_sequences[1].description == "DECOY_2 MASS=2077"
    assert out_sequences[1].seq == "ALEWLGFIRNKAOGYTT"


def test_random_sequences_invalid(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    output = "output.fasta"
    sequences = [SeqRecord(Seq("MFYHISLEHEILLHPRYFB"), id="P62487"),
                 SeqRecord(Seq("MEPDGTYEPGFVGIRFCQECNNML"), id="P36954")]
    aa_count = {'A': 17, 'C': 11, 'D': 21, 'E': 33}
    RandomSequences.parse_fasta = MagicMock(return_value=sequences)
    RandomSequences.aa_count = MagicMock(return_value=aa_count)
    RandomSequences.generate_sequence = MagicMock(
        side_effect=[Seq("EVKLVESGGGLVQPGGSL"), Seq("ALEWLGFIRNKAOGYTT")])
    with open(output, "w") as output_file:
        RandomSequences.random_sequences(fasta, output_file, invalid_sequence=True)
    RandomSequences.parse_fasta.assert_called_once_with(fasta)
    RandomSequences.aa_count.assert_called_once_with(sequences[1:])
    RandomSequences.generate_sequence.assert_called()
    assert RandomSequences.generate_sequence.call_count == 1
    assert os.path.isfile(output)
    out_sequences = []
    for record in SeqIO.parse(output, "fasta"):
        out_sequences.append(record)
    assert len(out_sequences) == 1
    assert out_sequences[0].name == "DECOY_1"
    assert out_sequences[0].description == "DECOY_1 MASS=1725"
    assert out_sequences[0].seq == "EVKLVESGGGLVQPGGSL"


def test_aa_count(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        sequences.append(record)
    aa_count = RandomSequences.aa_count(sequences)
    assert len(aa_count) == 20
    assert aa_count["A"] == 17
    assert aa_count["C"] == 11
    assert aa_count["D"] == 21
    assert aa_count["E"] == 33
    assert aa_count["F"] == 10
    assert aa_count["G"] == 14
    assert aa_count["H"] == 10
    assert aa_count["I"] == 20
    assert aa_count["K"] == 18
    assert aa_count["L"] == 25
    assert aa_count["M"] == 10
    assert aa_count["N"] == 10
    assert aa_count["P"] == 18
    assert aa_count["Q"] == 24
    assert aa_count["R"] == 23
    assert aa_count["S"] == 9
    assert aa_count["T"] == 22
    assert aa_count["V"] == 23
    assert aa_count["W"] == 2
    assert aa_count["Y"] == 15


def test_parse_fasta(testdir, mock_testclass):
    fasta = Path(__file__).parent.joinpath("P19388__P36954.fasta")
    sequences = RandomSequences.parse_fasta(fasta)
    assert len(sequences) == 2
    assert sequences[0].name == "sp|P19388|RPAB1_HUMAN"
    assert sequences[0].description == (
        "sp|P19388|RPAB1_HUMAN DNA-directed RNA polymerases I, II, and III subunit RPABC1 "
        "OS=Homo sapiens OX=9606 GN=POLR2E PE=1 SV=4")
    assert sequences[0].seq == ("MDDEEETYRLWKIRKTIMQLCHDRGYLVTQDELDQTLEEFKAQSGDKPSEGRPRRTDLTV"
                                "LVAHNDDPTDQMFVFFPEEPKVGIKTIKVYCQRMQEENITRALIVVQQGMTPSAKQSLVD"
                                "MAPKYILEQFLQQELLINITEHELVPEHVVMTKEEVTELLARYKLRENQLPRIQAGDPVA"
                                "RYFGIKRGQVVKIIRPSETAGRYITYRLVQ")
    assert sequences[1].name == "sp|P36954|RPB9_HUMAN"
    assert sequences[1].description == (
        "sp|P36954|RPB9_HUMAN DNA-directed RNA polymerase II subunit RPB9 "
        "OS=Homo sapiens OX=9606 GN=POLR2I PE=1 SV=1")
    assert sequences[1].seq == ("MEPDGTYEPGFVGIRFCQECNNMLYPKEDKENRILLYACRNCDYQQEADNSCIYVNKITH"
                                "EVDELTQIIADVSQDPTLPRTEDHPCQKCGHKEAVFFQSHSARAEDAMRLYYVCTAPHCG"
                                "HRWTE")


def test_generate_sequence(testdir, mock_testclass):
    aa_probabilities = {'A': 0.07011778585762114,
                        'C': 0.02299290897184481,
                        'D': 0.04738778390886124,
                        'E': 0.07108926900435926,
                        'F': 0.03647885419940203,
                        'G': 0.06571437851237638,
                        'H': 0.026215385371660117,
                        'I': 0.04337430410837212,
                        'K': 0.057355954163762546,
                        'L': 0.09962192302222815,
                        'M': 0.02132645590799044,
                        'N': 0.03586947171050206,
                        'P': 0.063154129352011,
                        'Q': 0.04769642534230582,
                        'R': 0.05637165484828641,
                        'S': 0.08330852029431542,
                        'T': 0.0535306789233716,
                        'V': 0.05962687392576805,
                        'W': 0.012147182286299164,
                        'Y': 0.02661690013746658}
    sequence = RandomSequences.generate_sequence(100000, aa_probabilities)
    assert len(sequence) == 100000
    seq_aa_counts = {}
    for aa in sequence:
        seq_aa_counts[aa] = seq_aa_counts[aa] + 1 if aa in seq_aa_counts else 1
    seq_aa_sum = sum(seq_aa_counts.values())
    seq_aa_probabilities = {aa: seq_aa_counts[aa] / seq_aa_sum for aa in seq_aa_counts}
    for aa in aa_probabilities:
        assert aa in seq_aa_probabilities
        assert abs(seq_aa_probabilities[aa] - aa_probabilities[
            aa]) < 0.005, f"{aa}: probability {seq_aa_probabilities[aa]} differs from {aa_probabilities[aa]}"
