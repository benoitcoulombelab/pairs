import os
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY

import pytest

from afpairs import ConsensusInterface
from afpairs.ConsensusInterface import ResiduePair, Residue, ConsensusResiduePair


@pytest.fixture
def mock_testclass():
    _consensus_interface = ConsensusInterface.consensus_interface
    _consensus_residue_pairs = ConsensusInterface.consensus_residue_pairs
    _parse_residue_pairs = ConsensusInterface.parse_residue_pairs
    _parse_alignment = ConsensusInterface.parse_alignment
    yield
    ConsensusInterface.consensus_interface = _consensus_interface
    ConsensusInterface.consensus_residue_pairs = _consensus_residue_pairs
    ConsensusInterface.parse_residue_pairs = _parse_residue_pairs
    ConsensusInterface.parse_alignment = _parse_alignment


def test_main(testdir, mock_testclass):
    residue_pairs_1 = "residue_pairs_1.txt"
    residue_pairs_2 = "residue_pairs_2.txt"
    open(residue_pairs_1, 'w').close()
    open(residue_pairs_2, 'w').close()
    ConsensusInterface.consensus_interface = MagicMock()
    ConsensusInterface.main(["-r", residue_pairs_1, residue_pairs_2])
    ConsensusInterface.consensus_interface.assert_called_once_with(
        residue_pair_files=ANY, output_file=ANY, name=r"([\w-]+)__([\w-]+)", consensus_ratio=0.5,
        baits_file=None, targets_file=None, alignment_format="clustal")
    residue_pair_in = ConsensusInterface.consensus_interface.call_args.kwargs["residue_pair_files"]
    assert len(residue_pair_in) == 2
    assert residue_pair_in[0].name == residue_pairs_1
    assert residue_pair_in[0].mode == "r"
    assert residue_pair_in[1].name == residue_pairs_2
    assert residue_pair_in[1].mode == "r"
    output_out = ConsensusInterface.consensus_interface.call_args.kwargs["output_file"]
    assert isinstance(output_out, TextIOWrapper)
    assert output_out.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
    residue_pairs_1 = "residue_pairs_1.txt"
    residue_pairs_2 = "residue_pairs_2.txt"
    open(residue_pairs_1, 'w').close()
    open(residue_pairs_2, 'w').close()
    output_file = "output.txt"
    name = r"(\w+)__(\w+)"
    consensus_ratio = 0.3
    baits_alignment = "baits_alignment.txt"
    open(baits_alignment, 'w').close()
    targets_alignment = "targets_alignment.txt"
    open(targets_alignment, 'w').close()
    alignment_format = "stockholm"
    ConsensusInterface.consensus_interface = MagicMock()
    ConsensusInterface.main(
        ["-o", output_file, "-n", name, "-c", str(consensus_ratio),
         "-b", baits_alignment, "-t", targets_alignment, "-f", alignment_format,
         "-r", residue_pairs_1, residue_pairs_2])
    ConsensusInterface.consensus_interface.assert_called_once_with(
        residue_pair_files=ANY, output_file=ANY, name=name, consensus_ratio=consensus_ratio,
        baits_file=ANY, targets_file=ANY, alignment_format=alignment_format)
    residue_pair_in = ConsensusInterface.consensus_interface.call_args.kwargs["residue_pair_files"]
    assert len(residue_pair_in) == 2
    assert residue_pair_in[0].name == residue_pairs_1
    assert residue_pair_in[0].mode == "r"
    assert residue_pair_in[1].name == residue_pairs_2
    assert residue_pair_in[1].mode == "r"
    output_out = ConsensusInterface.consensus_interface.call_args.kwargs["output_file"]
    assert output_out.name == output_file
    assert output_out.mode == "w"
    baits_in = ConsensusInterface.consensus_interface.call_args.kwargs["baits_file"]
    assert baits_in.name == baits_alignment
    assert baits_in.mode == "r"
    targets_in = ConsensusInterface.consensus_interface.call_args.kwargs["targets_file"]
    assert targets_in.name == targets_alignment
    assert targets_in.mode == "r"


def test_main_long_parameters(testdir, mock_testclass):
    residue_pairs_1 = "residue_pairs_1.txt"
    residue_pairs_2 = "residue_pairs_2.txt"
    open(residue_pairs_1, 'w').close()
    open(residue_pairs_2, 'w').close()
    output_file = "output.txt"
    name = r"(\w+)__(\w+)"
    consensus_ratio = 0.3
    baits_alignment = "baits_alignment.txt"
    open(baits_alignment, 'w').close()
    targets_alignment = "targets_alignment.txt"
    open(targets_alignment, 'w').close()
    alignment_format = "stockholm"
    ConsensusInterface.consensus_interface = MagicMock()
    ConsensusInterface.main(
        ["--output", output_file, "--name", name, "--consensus", str(consensus_ratio),
         "--baits", baits_alignment, "--targets", targets_alignment, "--format", alignment_format,
         "--residues", residue_pairs_1, residue_pairs_2])
    ConsensusInterface.consensus_interface.assert_called_once_with(
        residue_pair_files=ANY, output_file=ANY, name=name, consensus_ratio=consensus_ratio,
        baits_file=ANY, targets_file=ANY, alignment_format=alignment_format)
    residue_pair_in = ConsensusInterface.consensus_interface.call_args.kwargs["residue_pair_files"]
    assert len(residue_pair_in) == 2
    assert residue_pair_in[0].name == residue_pairs_1
    assert residue_pair_in[0].mode == "r"
    assert residue_pair_in[1].name == residue_pairs_2
    assert residue_pair_in[1].mode == "r"
    output_out = ConsensusInterface.consensus_interface.call_args.kwargs["output_file"]
    assert output_out.name == output_file
    assert output_out.mode == "w"
    baits_in = ConsensusInterface.consensus_interface.call_args.kwargs["baits_file"]
    assert baits_in.name == baits_alignment
    assert baits_in.mode == "r"
    targets_in = ConsensusInterface.consensus_interface.call_args.kwargs["targets_file"]
    assert targets_in.name == targets_alignment
    assert targets_in.mode == "r"


def test_consensus_interface(testdir, mock_testclass):
    residue_pairs_1_file = "HA_G4__HEMA_I78A3.txt"
    residue_pairs_2_file = "HEMA_I02A3__HEMA_I03A0.txt"
    open(residue_pairs_1_file, 'w').close()
    open(residue_pairs_2_file, 'w').close()
    output_file = "output.txt"
    baits_alignment_file = "baits.txt"
    open(baits_alignment_file, 'w').close()
    targets_alignment_file = "targets.txt"
    open(targets_alignment_file, 'w').close()
    baits = ["HA_G4", "HEMA_I02A3", "POLR2A"]
    targets = ["HEMA_I78A3", "HEMA_I03A0", "POLR2B"]
    baits_alignment = {bait: None for bait in baits}
    targets_alignment = {target: None for target in targets}
    residue_pairs_1 = [ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)]
    residue_pairs_2 = [ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)]
    consensus_residue_pairs = [
        ConsensusResiduePair(8, 131).append(
            "HA_G4", "HEMA_I78A3", ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)).append(
            "HEMA_I02A3", "HEMA_I78A3", ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)),
        ConsensusResiduePair(12, 275).append(
            "HA_G4", "HEMA_I78A3",
            ResiduePair(Residue("A", 11, "VAL"), Residue("B", 275, "ILE"), "Hydrophobic")).append(
            "HEMA_I02A3", "HEMA_I03A0", ResiduePair(Residue("A", 10, "PHE"), Residue("B", 274, "LEU"), None)),
        ConsensusResiduePair(19, 421).append(
            "HA_G4", "HEMA_I78A3",
            ResiduePair(Residue("A", 19, "THR"), Residue("B", 420, "THR"), "Charged")).append(
            "HEMA_I02A3", "HEMA_I03A0",
            ResiduePair(Residue("A", 18, "GLN"), Residue("B", 421, "ASP"), "Hydrogen")).append(
            "HA_G4", "HEMA_I03A0", ResiduePair(Residue("A", 19, "THR"), Residue("B", 420, "THR"), "Hydrogen"))]
    ConsensusInterface.parse_residue_pairs = MagicMock(side_effect=[residue_pairs_1, residue_pairs_2])
    ConsensusInterface.parse_alignment = MagicMock(side_effect=[baits_alignment, targets_alignment])
    ConsensusInterface.consensus_residue_pairs = MagicMock(return_value=consensus_residue_pairs)
    with open(residue_pairs_1_file, 'r') as residue_pairs_1_in, open(residue_pairs_2_file, 'r') as residue_pairs_2_in, \
            open(baits_alignment_file, 'r') as baits_in, open(targets_alignment_file, 'r') as targets_in, \
            open(output_file, 'w') as output_out:
        ConsensusInterface.consensus_interface(residue_pair_files=[residue_pairs_1_in, residue_pairs_2_in],
                                               output_file=output_out,
                                               baits_file=baits_in, targets_file=targets_in)
        ConsensusInterface.parse_alignment.assert_any_call(baits_in, "clustal")
        ConsensusInterface.parse_alignment.assert_any_call(targets_in, "clustal")
        ConsensusInterface.parse_residue_pairs.assert_any_call(residue_pairs_1_in)
        ConsensusInterface.parse_residue_pairs.assert_any_call(residue_pairs_2_in)
    ConsensusInterface.consensus_residue_pairs.assert_called_once_with(ANY, baits_alignment, targets_alignment)
    residue_pairs = ConsensusInterface.consensus_residue_pairs.call_args.args[0]
    assert len(residue_pairs) == 2
    assert ("HA_G4", "HEMA_I78A3") in residue_pairs
    assert residue_pairs[("HA_G4", "HEMA_I78A3")] == residue_pairs_1
    assert ("HEMA_I02A3", "HEMA_I03A0") in residue_pairs
    assert residue_pairs[("HEMA_I02A3", "HEMA_I03A0")] == residue_pairs_2
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "Bait residue index\tTarget residue index\tConsensus count\t" \
                                       "Baits\tTargets\t" \
                                       "Chain A\tResidue number A\tResidue name A\t" \
                                       "Chain B\tResidue number B\tResidue name B\tBond type (guess)\n"
        assert output_in.readline() == "8\t131\t2\t" \
                                       "HA_G4,HEMA_I02A3\tHEMA_I78A3,HEMA_I78A3\t" \
                                       "A,A\t7,6\tLEU,LEU\t" \
                                       "B,B\t131,131\tGLU,GLU\t,\n"
        assert output_in.readline() == "12\t275\t2\t" \
                                       "HA_G4,HEMA_I02A3\tHEMA_I78A3,HEMA_I03A0\t" \
                                       "A,A\t11,10\tVAL,PHE\t" \
                                       "B,B\t275,274\tILE,LEU\tHydrophobic,\n"
        assert output_in.readline() == "19\t421\t3\t" \
                                       "HA_G4,HEMA_I02A3,HA_G4\tHEMA_I78A3,HEMA_I03A0,HEMA_I03A0\t" \
                                       "A,A,A\t19,18,19\tTHR,GLN,THR\t" \
                                       "B,B,B\t420,421,420\tTHR,ASP,THR\tCharged,Hydrogen,Hydrogen\n"


def test_consensus_interface_no_alignment(testdir, mock_testclass):
    residue_pairs_1_file = "HA_G4__HEMA_I78A3.txt"
    residue_pairs_2_file = "HEMA_I02A3__HEMA_I03A0.txt"
    open(residue_pairs_1_file, 'w').close()
    open(residue_pairs_2_file, 'w').close()
    output_file = "output.txt"
    residue_pairs_1 = [ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)]
    residue_pairs_2 = [ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)]
    consensus_residue_pairs = [
        ConsensusResiduePair(8, 131).append(
            "HA_G4", "HEMA_I78A3", ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)).append(
            "HEMA_I02A3", "HEMA_I78A3", ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)),
        ConsensusResiduePair(12, 275).append(
            "HA_G4", "HEMA_I78A3",
            ResiduePair(Residue("A", 11, "VAL"), Residue("B", 275, "ILE"), "Hydrophobic")).append(
            "HEMA_I02A3", "HEMA_I03A0", ResiduePair(Residue("A", 10, "PHE"), Residue("B", 274, "LEU"), None)),
        ConsensusResiduePair(19, 421).append(
            "HA_G4", "HEMA_I78A3",
            ResiduePair(Residue("A", 19, "THR"), Residue("B", 420, "THR"), "Charged")).append(
            "HEMA_I02A3", "HEMA_I03A0",
            ResiduePair(Residue("A", 18, "GLN"), Residue("B", 421, "ASP"), "Hydrogen")).append(
            "HA_G4", "HEMA_I03A0", ResiduePair(Residue("A", 19, "THR"), Residue("B", 420, "THR"), "Hydrogen"))]
    ConsensusInterface.parse_residue_pairs = MagicMock(side_effect=[residue_pairs_1, residue_pairs_2])
    ConsensusInterface.parse_alignment = MagicMock()
    ConsensusInterface.consensus_residue_pairs = MagicMock(return_value=consensus_residue_pairs)
    with open(residue_pairs_1_file, 'r') as residue_pairs_1_in, open(residue_pairs_2_file, 'r') as residue_pairs_2_in, \
            open(output_file, 'w') as output_out:
        ConsensusInterface.consensus_interface(residue_pair_files=[residue_pairs_1_in, residue_pairs_2_in],
                                               output_file=output_out)
        ConsensusInterface.parse_alignment.assert_not_called()
        ConsensusInterface.parse_residue_pairs.assert_any_call(residue_pairs_1_in)
        ConsensusInterface.parse_residue_pairs.assert_any_call(residue_pairs_2_in)
    ConsensusInterface.consensus_residue_pairs.assert_called_once_with(ANY, None, None)
    residue_pairs = ConsensusInterface.consensus_residue_pairs.call_args.args[0]
    assert len(residue_pairs) == 2
    assert ("HA_G4", "HEMA_I78A3") in residue_pairs
    assert residue_pairs[("HA_G4", "HEMA_I78A3")] == residue_pairs_1
    assert ("HEMA_I02A3", "HEMA_I03A0") in residue_pairs
    assert residue_pairs[("HEMA_I02A3", "HEMA_I03A0")] == residue_pairs_2
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "Bait residue index\tTarget residue index\tConsensus count\t" \
                                       "Baits\tTargets\t" \
                                       "Chain A\tResidue number A\tResidue name A\t" \
                                       "Chain B\tResidue number B\tResidue name B\tBond type (guess)\n"
        assert output_in.readline() == "8\t131\t2\t" \
                                       "HA_G4,HEMA_I02A3\tHEMA_I78A3,HEMA_I78A3\t" \
                                       "A,A\t7,6\tLEU,LEU\t" \
                                       "B,B\t131,131\tGLU,GLU\t,\n"
        assert output_in.readline() == "12\t275\t2\t" \
                                       "HA_G4,HEMA_I02A3\tHEMA_I78A3,HEMA_I03A0\t" \
                                       "A,A\t11,10\tVAL,PHE\t" \
                                       "B,B\t275,274\tILE,LEU\tHydrophobic,\n"
        assert output_in.readline() == "19\t421\t3\t" \
                                       "HA_G4,HEMA_I02A3,HA_G4\tHEMA_I78A3,HEMA_I03A0,HEMA_I03A0\t" \
                                       "A,A,A\t19,18,19\tTHR,GLN,THR\t" \
                                       "B,B,B\t420,421,420\tTHR,ASP,THR\tCharged,Hydrogen,Hydrogen\n"


def test_consensus_interface_missing_bait_in_alignment(testdir, mock_testclass):
    residue_pairs_1_file = "HA_G4__HEMA_I78A3.txt"
    residue_pairs_2_file = "HEMA_I02A3__HEMA_I03A0.txt"
    open(residue_pairs_1_file, 'w').close()
    open(residue_pairs_2_file, 'w').close()
    output_file = "output.txt"
    baits_alignment_file = "baits.txt"
    open(baits_alignment_file, 'w').close()
    targets_alignment_file = "targets.txt"
    open(targets_alignment_file, 'w').close()
    baits = ["HA_G4", "POLR2A"]
    targets = ["HEMA_I78A3", "HEMA_I03A0", "POLR2B"]
    baits_alignment = {bait: None for bait in baits}
    targets_alignment = {target: None for target in targets}
    residue_pairs_1 = [ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)]
    residue_pairs_2 = [ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)]
    ConsensusInterface.parse_residue_pairs = MagicMock(side_effect=[residue_pairs_1, residue_pairs_2])
    ConsensusInterface.parse_alignment = MagicMock(side_effect=[baits_alignment, targets_alignment])
    ConsensusInterface.consensus_residue_pairs = MagicMock()
    with open(residue_pairs_1_file, 'r') as residue_pairs_1_in, open(residue_pairs_2_file, 'r') as residue_pairs_2_in, \
            open(baits_alignment_file, 'r') as baits_in, open(targets_alignment_file, 'r') as targets_in, \
            open(output_file, 'w') as output_out, pytest.raises(AssertionError) as e_info:
        ConsensusInterface.consensus_interface(residue_pair_files=[residue_pairs_1_in, residue_pairs_2_in],
                                               output_file=output_out,
                                               baits_file=baits_in, targets_file=targets_in)
    ConsensusInterface.consensus_residue_pairs.assert_not_called()
    assert os.path.getsize(output_file) == 0


def test_consensus_interface_missing_target_in_alignment(testdir, mock_testclass):
    residue_pairs_1_file = "HA_G4__HEMA_I78A3.txt"
    residue_pairs_2_file = "HEMA_I02A3__HEMA_I03A0.txt"
    open(residue_pairs_1_file, 'w').close()
    open(residue_pairs_2_file, 'w').close()
    output_file = "output.txt"
    baits_alignment_file = "baits.txt"
    open(baits_alignment_file, 'w').close()
    targets_alignment_file = "targets.txt"
    open(targets_alignment_file, 'w').close()
    baits = ["HA_G4", "HEMA_I02A3", "POLR2A"]
    targets = ["HEMA_I78A3", "POLR2B"]
    baits_alignment = {bait: None for bait in baits}
    targets_alignment = {target: None for target in targets}
    residue_pairs_1 = [ResiduePair(Residue("A", 7, "LEU"), Residue("B", 131, "GLU"), None)]
    residue_pairs_2 = [ResiduePair(Residue("A", 6, "LEU"), Residue("B", 131, "GLU"), None)]
    ConsensusInterface.parse_residue_pairs = MagicMock(side_effect=[residue_pairs_1, residue_pairs_2])
    ConsensusInterface.parse_alignment = MagicMock(side_effect=[baits_alignment, targets_alignment])
    ConsensusInterface.consensus_residue_pairs = MagicMock()
    with open(residue_pairs_1_file, 'r') as residue_pairs_1_in, open(residue_pairs_2_file, 'r') as residue_pairs_2_in, \
            open(baits_alignment_file, 'r') as baits_in, open(targets_alignment_file, 'r') as targets_in, \
            open(output_file, 'w') as output_out, pytest.raises(AssertionError) as e_info:
        ConsensusInterface.consensus_interface(residue_pair_files=[residue_pairs_1_in, residue_pairs_2_in],
                                               output_file=output_out,
                                               baits_file=baits_in, targets_file=targets_in)
    ConsensusInterface.consensus_residue_pairs.assert_not_called()
    assert os.path.getsize(output_file) == 0


def test_consensus_residue_pairs(testdir, mock_testclass):
    baits = ["HA_G4", "HEMA_I78A3", "HEMA_I02A3"]
    target = "POLR2A"
    residue_pairs = {(bait, target): [] for bait in baits}
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 3, "ILE"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 5, "PHE"), Residue("B", 280, "TYR"), "Hydrophobic"))
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 534, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[1], target)].append(
        ResiduePair(Residue("A", 3, "MET"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[1], target)].append(
        ResiduePair(Residue("A", 534, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[2], target)].append(
        ResiduePair(Residue("A", 3, "LYS"), Residue("B", 130, "GLU"), "Hydrogen"))
    residue_pairs[(baits[2], target)].append(
        ResiduePair(Residue("A", 8, "ILE"), Residue("B", 131, "GLU"), None))
    consensus_residue_pairs = ConsensusInterface.consensus_residue_pairs(residue_pairs)
    assert len(consensus_residue_pairs) == 4
    consensus_residue_pair = consensus_residue_pairs[0]
    assert consensus_residue_pair.bait_residue_index == 3
    assert consensus_residue_pair.target_residue_index == 130
    assert len(consensus_residue_pair.residue_pairs) == 3
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][0]
    assert (baits[1], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], target)] == residue_pairs[(baits[1], target)][0]
    assert (baits[2], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[2], target)] == residue_pairs[(baits[2], target)][0]
    consensus_residue_pair = consensus_residue_pairs[1]
    assert consensus_residue_pair.bait_residue_index == 5
    assert consensus_residue_pair.target_residue_index == 280
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][1]
    consensus_residue_pair = consensus_residue_pairs[2]
    assert consensus_residue_pair.bait_residue_index == 534
    assert consensus_residue_pair.target_residue_index == 455
    assert len(consensus_residue_pair.residue_pairs) == 2
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][2]
    assert (baits[1], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], target)] == residue_pairs[(baits[1], target)][1]
    consensus_residue_pair = consensus_residue_pairs[3]
    assert consensus_residue_pair.bait_residue_index == 8
    assert consensus_residue_pair.target_residue_index == 131
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[2], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[2], target)] == residue_pairs[(baits[2], target)][1]


def test_consensus_residue_pairs_bait_align(testdir, mock_testclass):
    baits = ["HA_G4", "HEMA_I78A3", "HEMA_I02A3"]
    target = "POLR2A"
    residue_pairs = {(bait, target): [] for bait in baits}
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 3, "ILE"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 5, "PHE"), Residue("B", 280, "TYR"), "Hydrophobic"))
    residue_pairs[(baits[0], target)].append(
        ResiduePair(Residue("A", 530, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[1], target)].append(
        ResiduePair(Residue("A", 1, "MET"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[1], target)].append(
        ResiduePair(Residue("A", 526, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[2], target)].append(
        ResiduePair(Residue("A", 2, "LYS"), Residue("B", 130, "GLU"), "Hydrogen"))
    residue_pairs[(baits[2], target)].append(
        ResiduePair(Residue("A", 7, "ILE"), Residue("B", 131, "GLU"), None))
    alignment_file = Path(__file__).parent.joinpath("clustal_alignment.txt")
    with open(alignment_file, "r") as alignment_in:
        alignment = ConsensusInterface.parse_alignment(alignment_in, "clustal")
    consensus_residue_pairs = ConsensusInterface.consensus_residue_pairs(residue_pairs, baits=alignment)
    assert len(consensus_residue_pairs) == 4
    consensus_residue_pair = consensus_residue_pairs[0]
    assert consensus_residue_pair.bait_residue_index == 3
    assert consensus_residue_pair.target_residue_index == 130
    assert len(consensus_residue_pair.residue_pairs) == 3
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][0]
    assert (baits[1], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], target)] == residue_pairs[(baits[1], target)][0]
    assert (baits[2], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[2], target)] == residue_pairs[(baits[2], target)][0]
    consensus_residue_pair = consensus_residue_pairs[1]
    assert consensus_residue_pair.bait_residue_index == 5
    assert consensus_residue_pair.target_residue_index == 280
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][1]
    consensus_residue_pair = consensus_residue_pairs[2]
    assert consensus_residue_pair.bait_residue_index == 534
    assert consensus_residue_pair.target_residue_index == 455
    assert len(consensus_residue_pair.residue_pairs) == 2
    assert (baits[0], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], target)] == residue_pairs[(baits[0], target)][2]
    assert (baits[1], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], target)] == residue_pairs[(baits[1], target)][1]
    consensus_residue_pair = consensus_residue_pairs[3]
    assert consensus_residue_pair.bait_residue_index == 8
    assert consensus_residue_pair.target_residue_index == 131
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[2], target) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[2], target)] == residue_pairs[(baits[2], target)][1]


def test_consensus_residue_pairs_target_align(testdir, mock_testclass):
    bait = "POLR2A"
    targets = ["HA_G4", "HEMA_I78A3", "HEMA_I02A3"]
    residue_pairs = {(bait, target): [] for target in targets}
    residue_pairs[(bait, targets[0])].append(
        ResiduePair(Residue("B", 130, "GLU"), Residue("A", 3, "ILE"), None))
    residue_pairs[(bait, targets[0])].append(
        ResiduePair(Residue("B", 280, "TYR"), Residue("A", 5, "PHE"), "Hydrophobic"))
    residue_pairs[(bait, targets[0])].append(
        ResiduePair(Residue("B", 455, "LYS"), Residue("A", 530, "TYR"), "Hydrogen"))
    residue_pairs[(bait, targets[1])].append(
        ResiduePair(Residue("B", 130, "GLU"), Residue("A", 1, "MET"), None))
    residue_pairs[(bait, targets[1])].append(
        ResiduePair(Residue("B", 455, "LYS"), Residue("A", 526, "TYR"), "Hydrogen"))
    residue_pairs[(bait, targets[2])].append(
        ResiduePair(Residue("B", 130, "GLU"), Residue("A", 2, "LYS"), "Hydrogen"))
    residue_pairs[(bait, targets[2])].append(
        ResiduePair(Residue("B", 131, "GLU"), Residue("A", 7, "ILE"), None))
    alignment_file = Path(__file__).parent.joinpath("clustal_alignment.txt")
    with open(alignment_file, "r") as alignment_in:
        alignment = ConsensusInterface.parse_alignment(alignment_in, "clustal")
    consensus_residue_pairs = ConsensusInterface.consensus_residue_pairs(residue_pairs, targets=alignment)
    assert len(consensus_residue_pairs) == 4
    consensus_residue_pair = consensus_residue_pairs[0]
    assert consensus_residue_pair.bait_residue_index == 130
    assert consensus_residue_pair.target_residue_index == 3
    assert len(consensus_residue_pair.residue_pairs) == 3
    assert (bait, targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[0])] == residue_pairs[(bait, targets[0])][0]
    assert (bait, targets[1]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[1])] == residue_pairs[(bait, targets[1])][0]
    assert (bait, targets[2]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[2])] == residue_pairs[(bait, targets[2])][0]
    consensus_residue_pair = consensus_residue_pairs[1]
    assert consensus_residue_pair.bait_residue_index == 280
    assert consensus_residue_pair.target_residue_index == 5
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (bait, targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[0])] == residue_pairs[(bait, targets[0])][1]
    consensus_residue_pair = consensus_residue_pairs[2]
    assert consensus_residue_pair.bait_residue_index == 455
    assert consensus_residue_pair.target_residue_index == 534
    assert len(consensus_residue_pair.residue_pairs) == 2
    assert (bait, targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[0])] == residue_pairs[(bait, targets[0])][2]
    assert (bait, targets[1]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[1])] == residue_pairs[(bait, targets[1])][1]
    consensus_residue_pair = consensus_residue_pairs[3]
    assert consensus_residue_pair.bait_residue_index == 131
    assert consensus_residue_pair.target_residue_index == 8
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (bait, targets[2]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(bait, targets[2])] == residue_pairs[(bait, targets[2])][1]


def test_consensus_residue_pairs_bait_align(testdir, mock_testclass):
    baits = ["HA_G4", "HEMA_I02A3"]
    targets = ["HEMA_I78A3", "HEMA_I03A0"]
    residue_pairs = {(bait, target): [] for bait in baits for target in targets}
    residue_pairs[(baits[0], targets[0])].append(
        ResiduePair(Residue("A", 3, "ILE"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[0], targets[0])].append(
        ResiduePair(Residue("A", 5, "PHE"), Residue("B", 280, "TYR"), "Hydrophobic"))
    residue_pairs[(baits[0], targets[0])].append(
        ResiduePair(Residue("A", 530, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[1], targets[0])].append(
        ResiduePair(Residue("A", 2, "MET"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[1], targets[0])].append(
        ResiduePair(Residue("A", 532, "TYR"), Residue("B", 455, "LYS"), "Hydrogen"))
    residue_pairs[(baits[0], targets[1])].append(
        ResiduePair(Residue("A", 3, "LYS"), Residue("B", 131, "GLU"), "Hydrogen"))
    residue_pairs[(baits[0], targets[1])].append(
        ResiduePair(Residue("A", 7, "ILE"), Residue("B", 130, "GLU"), None))
    residue_pairs[(baits[1], targets[1])].append(
        ResiduePair(Residue("A", 532, "LYS"), Residue("B", 461, "GLU"), "Hydrogen"))
    alignment_file = Path(__file__).parent.joinpath("clustal_alignment.txt")
    with open(alignment_file, "r") as alignment_in:
        alignment = ConsensusInterface.parse_alignment(alignment_in, "clustal")
    consensus_residue_pairs = ConsensusInterface.consensus_residue_pairs(
        residue_pairs, baits=alignment, targets=alignment)
    assert len(consensus_residue_pairs) == 4
    consensus_residue_pair = consensus_residue_pairs[0]
    assert consensus_residue_pair.bait_residue_index == 3
    assert consensus_residue_pair.target_residue_index == 132
    assert len(consensus_residue_pair.residue_pairs) == 3
    assert (baits[0], targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], targets[0])] == residue_pairs[(baits[0], targets[0])][0]
    assert (baits[1], targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], targets[0])] == residue_pairs[(baits[1], targets[0])][0]
    assert (baits[0], targets[1]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], targets[1])] == residue_pairs[(baits[0], targets[1])][0]
    consensus_residue_pair = consensus_residue_pairs[1]
    assert consensus_residue_pair.bait_residue_index == 5
    assert consensus_residue_pair.target_residue_index == 284
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[0], targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], targets[0])] == residue_pairs[(baits[0], targets[0])][1]
    consensus_residue_pair = consensus_residue_pairs[2]
    assert consensus_residue_pair.bait_residue_index == 534
    assert consensus_residue_pair.target_residue_index == 463
    assert len(consensus_residue_pair.residue_pairs) == 3
    assert (baits[0], targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], targets[0])] == residue_pairs[(baits[0], targets[0])][2]
    assert (baits[1], targets[0]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], targets[0])] == residue_pairs[(baits[1], targets[0])][1]
    assert (baits[1], targets[1]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[1], targets[1])] == residue_pairs[(baits[1], targets[1])][0]
    consensus_residue_pair = consensus_residue_pairs[3]
    assert consensus_residue_pair.bait_residue_index == 7
    assert consensus_residue_pair.target_residue_index == 131
    assert len(consensus_residue_pair.residue_pairs) == 1
    assert (baits[0], targets[1]) in consensus_residue_pair.residue_pairs
    assert consensus_residue_pair.residue_pairs[(baits[0], targets[1])] == residue_pairs[(baits[0], targets[1])][1]


def test_parse_residue_pairs(mock_testclass):
    residue_pairs_file = Path(__file__).parent.joinpath("residue_pairs.txt")
    with open(residue_pairs_file, "r") as residue_pairs_in:
        residue_pairs = ConsensusInterface.parse_residue_pairs(residue_pairs_in)
    assert len(residue_pairs) == 318
    residue_pair = residue_pairs[0]
    assert residue_pair.first.chain == "A"
    assert residue_pair.first.index == 3
    assert residue_pair.first.name == "ILE"
    assert residue_pair.second.chain == "B"
    assert residue_pair.second.index == 130
    assert residue_pair.second.name == "GLU"
    assert residue_pair.bond is None
    residue_pair = residue_pairs[4]
    assert residue_pair.first.chain == "A"
    assert residue_pair.first.index == 5
    assert residue_pair.first.name == "PHE"
    assert residue_pair.second.chain == "B"
    assert residue_pair.second.index == 280
    assert residue_pair.second.name == "TYR"
    assert residue_pair.bond == "Hydrophobic"
    residue_pair = residue_pairs[317]
    assert residue_pair.first.chain == "A"
    assert residue_pair.first.index == 530
    assert residue_pair.first.name == "TYR"
    assert residue_pair.second.chain == "B"
    assert residue_pair.second.index == 455
    assert residue_pair.second.name == "LYS"
    assert residue_pair.bond == "Hydrogen"


def test_parse_alignment(testdir, mock_testclass):
    alignment_file = Path(__file__).parent.joinpath("clustal_alignment.txt")
    with open(alignment_file, "r") as alignment_in:
        alignment = ConsensusInterface.parse_alignment(alignment_in, "clustal")
    assert "HA_G4" in alignment
    sequence = alignment["HA_G4"]
    assert sequence.position(10) == 10
    assert sequence.position(356) == 360
    assert "HEMA_I78A3" in alignment
    sequence = alignment["HEMA_I78A3"]
    assert sequence.position(10) == 12
    assert sequence.position(356) == 364
    assert "HEMA_I02A3" in alignment
    sequence = alignment["HEMA_I02A3"]
    assert sequence.position(10) == 11
    assert sequence.position(356) == 358
    assert "HEMA_I03A0" in alignment
    sequence = alignment["HEMA_I03A0"]
    assert sequence.position(10) == 11
    assert sequence.position(356) == 358
