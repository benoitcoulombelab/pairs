from io import TextIOWrapper
from unittest.mock import MagicMock, ANY

import pytest

from afpairs import IdConvert


@pytest.fixture
def mock_testclass():
    _id_convert = IdConvert.id_convert
    _parse_mapping = IdConvert.parse_mapping
    yield
    IdConvert.id_convert = _id_convert
    IdConvert.parse_mapping = _parse_mapping


def test_main(testdir, mock_testclass):
    mapping_file = "mapping.txt"
    open(mapping_file, 'w').close()
    IdConvert.id_convert = MagicMock()
    IdConvert.main()
    IdConvert.id_convert.assert_called_once_with(
        input_file=ANY, output_file=ANY, mapping_file=ANY, id_column=0,
        mapping_source_column=0, mapping_converted_column=1)
    assert IdConvert.id_convert.call_args.kwargs["input_file"].name == __file__
    assert IdConvert.id_convert.call_args.kwargs["input_file"].mode == "r"
    assert isinstance(IdConvert.id_convert.call_args.kwargs["output_file"], TextIOWrapper)
    assert IdConvert.id_convert.call_args.kwargs["output_file"].mode in ["r+", "w"]
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].name == mapping_file
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].mode == "r"


def test_main_parameters(testdir, mock_testclass):
    input_file = "input.txt"
    output_file = "output.txt"
    mapping_file = "mapping_file.txt"
    open(input_file, 'w').close()
    open(mapping_file, 'w').close()
    IdConvert.id_convert = MagicMock()
    IdConvert.main(["-i", "2", "-m", mapping_file, "-s", "3", "-c", "4", input_file, output_file])
    IdConvert.id_convert.assert_called_once_with(
        input_file=ANY, output_file=ANY, mapping_file=ANY, id_column=1,
        mapping_source_column=2, mapping_converted_column=3)
    assert IdConvert.id_convert.call_args.kwargs["input_file"].name == input_file
    assert IdConvert.id_convert.call_args.kwargs["input_file"].mode == "r"
    assert IdConvert.id_convert.call_args.kwargs["output_file"].name == output_file
    assert IdConvert.id_convert.call_args.kwargs["output_file"].mode == "w"
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].name == mapping_file
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].mode == "r"


def test_main_long_parameters(testdir, mock_testclass):
    input_file = "input.txt"
    output_file = "output.txt"
    mapping_file = "mapping_file.txt"
    open(input_file, 'w').close()
    open(mapping_file, 'w').close()
    IdConvert.id_convert = MagicMock()
    IdConvert.main(
        ["--id_column", "2", "--mapping", mapping_file, "--source_column", "3",
         "--converted_column", "4", input_file, output_file])
    IdConvert.id_convert.assert_called_once_with(
        input_file=ANY, output_file=ANY, mapping_file=ANY, id_column=1,
        mapping_source_column=2, mapping_converted_column=3)
    assert IdConvert.id_convert.call_args.kwargs["input_file"].name == input_file
    assert IdConvert.id_convert.call_args.kwargs["input_file"].mode == "r"
    assert IdConvert.id_convert.call_args.kwargs["output_file"].name == output_file
    assert IdConvert.id_convert.call_args.kwargs["output_file"].mode == "w"
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].name == mapping_file
    assert IdConvert.id_convert.call_args.kwargs["mapping_file"].mode == "r"


def test_id_convert(testdir, mock_testclass):
    input_file = "input.txt"
    output_file = "output.txt"
    mapping_file = "mapping.txt"
    mappings = {"RPB1_HUMAN": "POLR2A", "RPB2_HUMAN": "POLR2B"}
    IdConvert.parse_mapping = MagicMock(return_value=mappings)
    with open(input_file, 'w') as input_out:
        input_out.write("RPB1_HUMAN\t1\n")
        input_out.write("RPB2_HUMAN\t2\n")
    open(mapping_file, 'w').close()
    with open(input_file, 'r') as input_in, open(mapping_file, 'r') as mapping_in, \
            open(output_file, 'w') as output_out:
        IdConvert.id_convert(input_file=input_in, output_file=output_out, mapping_file=mapping_in)
    IdConvert.parse_mapping.assert_called_once_with(
        mapping_file=mapping_in, source_column=0, converted_column=1)
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "RPB1_HUMAN\tPOLR2A\t1\n"
        assert output_in.readline() == "RPB2_HUMAN\tPOLR2B\t2\n"


def test_id_convert_header(testdir, mock_testclass):
    input_file = "input.txt"
    output_file = "output.txt"
    mapping_file = "mapping.txt"
    mappings = {"RPB1_HUMAN": "POLR2A", "RPB2_HUMAN": "POLR2B"}
    IdConvert.parse_mapping = MagicMock(return_value=mappings)
    with open(input_file, 'w') as input_out:
        input_out.write("Protein\tScore\n")
        input_out.write("RPB1_HUMAN\t1\n")
        input_out.write("RPB2_HUMAN\t2\n")
    open(mapping_file, 'w').close()
    with open(input_file, 'r') as input_in, open(mapping_file, 'r') as mapping_in, \
            open(output_file, 'w') as output_out:
        IdConvert.id_convert(input_file=input_in, output_file=output_out, mapping_file=mapping_in)
    IdConvert.parse_mapping.assert_called_once_with(
        mapping_file=mapping_in, source_column=0, converted_column=1)
    with open(output_file, 'r') as output_in:
        assert output_in.readline() == "Protein\t\tScore\n"
        assert output_in.readline() == "RPB1_HUMAN\tPOLR2A\t1\n"
        assert output_in.readline() == "RPB2_HUMAN\tPOLR2B\t2\n"


def test_parse_mapping(testdir, mock_testclass):
    mapping_file = "mapping_file.txt"
    with open(mapping_file, 'w') as mapping_out:
        mapping_out.write("RPB1_HUMAN\tPOLR2A\n")
        mapping_out.write("RPB2_HUMAN\tPOLR2B\n")
    with open(mapping_file, 'r') as mapping_in:
        mappings = IdConvert.parse_mapping(mapping_file=mapping_in)
    assert "RPB1_HUMAN" in mappings
    assert mappings["RPB1_HUMAN"] == "POLR2A"
    assert "RPB2_HUMAN" in mappings
    assert mappings["RPB2_HUMAN"] == "POLR2B"
