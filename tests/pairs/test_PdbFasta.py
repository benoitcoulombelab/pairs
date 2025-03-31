from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, ANY, patch

import pytest

from pairs import PdbFasta


@pytest.fixture
def mock_testclass():
  _pdb_fasta = PdbFasta.pdb_fasta
  yield
  PdbFasta.pdb_fasta = _pdb_fasta


def test_main(testdir, mock_testclass):
  PdbFasta.pdb_fasta = MagicMock()
  stdin_file = "stdin.txt"
  open(stdin_file, 'w').close()
  with open(stdin_file, 'r') as stdin_in, patch('sys.stdin', stdin_in):
    PdbFasta.main()
  PdbFasta.pdb_fasta.assert_called_once_with(pdb=ANY, output=ANY, prefix="")
  pdb = PdbFasta.pdb_fasta.call_args.kwargs['pdb']
  assert isinstance(pdb, TextIOWrapper)
  assert pdb.mode == "r"
  output_file = PdbFasta.pdb_fasta.call_args.kwargs["output"]
  assert isinstance(output_file, TextIOWrapper)
  assert output_file.mode in ["r+", "w"]


def test_main_parameters(testdir, mock_testclass):
  pdb = "my.pdb"
  open(pdb, 'w').close()
  output = "output.txt"
  prefix = "my_prefix"
  PdbFasta.pdb_fasta = MagicMock()
  PdbFasta.main(["-p", prefix, "-o", output, pdb])
  PdbFasta.pdb_fasta.assert_called_once_with(pdb=ANY, output=ANY, prefix=prefix)
  pdb_input = PdbFasta.pdb_fasta.call_args.kwargs['pdb']
  assert isinstance(pdb_input, TextIOWrapper)
  assert pdb_input.mode == "r"
  assert pdb_input.name == pdb
  output_file = PdbFasta.pdb_fasta.call_args.kwargs["output"]
  assert isinstance(output_file, TextIOWrapper)
  assert output_file.mode == "w"
  assert output_file.name == output


def test_main_long_parameters(testdir, mock_testclass):
  pdb = "my.pdb"
  open(pdb, 'w').close()
  output = "output.txt"
  prefix = "my_prefix"
  PdbFasta.pdb_fasta = MagicMock()
  PdbFasta.main(["--prefix", prefix, "--output", output, pdb])
  PdbFasta.pdb_fasta.assert_called_once_with(pdb=ANY, output=ANY, prefix=prefix)
  pdb_input = PdbFasta.pdb_fasta.call_args.kwargs['pdb']
  assert isinstance(pdb_input, TextIOWrapper)
  assert pdb_input.mode == "r"
  assert pdb_input.name == pdb
  output_file = PdbFasta.pdb_fasta.call_args.kwargs["output"]
  assert isinstance(output_file, TextIOWrapper)
  assert output_file.mode == "w"
  assert output_file.name == output


def test_pdb_fasta(testdir, mock_testclass):
  pdb = Path(__file__).parent.joinpath("POLR2A_POLR2B_ranked_0.pdb")
  output = "output.txt"
  with open(pdb, "r") as pdb_in, open(output, "w") as output_out:
    PdbFasta.pdb_fasta(pdb_in, output_out)
  with open(output, "r") as output_in:
    assert output_in.readline() == ">A\n"
    assert output_in.readline() == "MHGGGPPSGDSACPLRTIKRVQFGVLSPDELKRMSVTEGGIKYPETTEGGRPKLGGLMDP\n"
    for i in range(0, 31):
      output_in.readline()
    assert output_in.readline() == "PKYSPTSPTYSPTSPKGSTYSPTSPGYSPTSPTYSLTSPAISPDDSDEEN\n"
    assert output_in.readline() == ">B\n"
    assert output_in.readline() == "MYDADEDMQYDEDDDEITPDLWQEACWIVISSYFDEKGLVRQQLDSFDEFIQMSVQRIVE\n"
    for i in range(0, 18):
      output_in.readline()
    assert output_in.readline() == "RNKTQISLVRMPYACKLLFQELMSMSIAPRMMSV\n"


def test_pdb_fasta_parameters(testdir, mock_testclass):
  pdb = Path(__file__).parent.joinpath("FAB_5_3__HVM13_MOUSE_ranked_0.pdb")
  output = "output.txt"
  prefix = "FAB_5_3__"
  with open(pdb, "r") as pdb_in, open(output, "w") as output_out:
    PdbFasta.pdb_fasta(pdb_in, output_out, prefix=prefix)
  with open(output, "r") as output_in:
    assert output_in.readline() == f">{prefix}A\n"
    assert output_in.readline() == "DVVMTQSPISLPVTPGEPASISCRSSQSLLFSNGYNYLDWYLQKPGQSPQLLIYLGSNRA\n"
    assert output_in.readline() == "SGVPDRFSGSGSGTDFTLQISRVEAEDVGVYYCMQALQTPLTFGGGTKVEIKRTVAAPSV\n"
    assert output_in.readline() == "FIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSL\n"
    assert output_in.readline() == "SSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGECGAAHHHHHHGAADYKDDDDKG\n"
    assert output_in.readline() == "AA\n"
    assert output_in.readline() == f">{prefix}B\n"
    assert output_in.readline() == "EVQLVQSGGGLVKPGGSLRLSCGTSGFSLRDSHMSWIRQAPGKGLEWIAYISRSGRNTIY\n"
    assert output_in.readline() == "PDSAKGRFTISRDNAKNTVSLQINSLRVEDTAVYYCAREKGDYYESGGAFDVWGQGTKVT\n"
    assert output_in.readline() == "VSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVL\n"
    assert output_in.readline() == "QSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDAAAEQKLISEEDLNAAA\n"
    assert output_in.readline() == f">{prefix}C\n"
    assert output_in.readline() == "EVQLQQSGPELVKPGASVKMSCKASGYTFTDYYMKWVKQSHGKSLEWIGDINPNNGGTSY\n"
    assert output_in.readline() == "NQKFKGKATLTVDKSSSTAYMQLNSLTSEDSAVYYCARDRYWYFDVWGAGTTVTVSS\n"
