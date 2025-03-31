import argparse
import sys
from typing import TextIO


def main(argv: list[str] = None):
  parser = argparse.ArgumentParser(description="Converts identifiers.")
  parser.add_argument('input', nargs='?', type=argparse.FileType('r'),
                      default=sys.stdin,
                      help="Tab delimited text file containing id to convert")
  parser.add_argument('-i', '--id_column', type=int, default='1',
                      help="Column index of id in input file - 1 means first column  (default: %(default)s)")
  parser.add_argument('-m', '--mapping', type=argparse.FileType('r'),
                      default="mapping.txt",
                      help="Tab delimited text file containing source ids and converted ids  (default: %(default)s)")
  parser.add_argument('-s', '--source_column', type=int, default='1',
                      help="Column index of source ids in mapping file - 1 means first column of file" +
                           "   (default: %(default)s)")
  parser.add_argument('-c', '--converted_column', type=int, default='2',
                      help="Column index of converted ids in mapping file - 1 means first column of file" +
                           "   (default: %(default)s)")
  parser.add_argument('output', nargs='?', type=argparse.FileType('w'),
                      default=sys.stdout,
                      help="Tab delimited output file")

  args = parser.parse_args(argv)
  id_convert(input_file=args.input, output_file=args.output,
             mapping_file=args.mapping,
             id_column=args.id_column - 1,
             mapping_source_column=args.source_column - 1,
             mapping_converted_column=args.converted_column - 1)


def id_convert(input_file: TextIO, output_file: TextIO, mapping_file: TextIO,
    id_column: int = 0,
    mapping_source_column: int = 0, mapping_converted_column: int = 1):
  """
  Converts id to another id.

  :param input_file: text delimited input file
  :param output_file: text delimited input file
  :param mapping_file: text delimited file containing mapping from source id to converted id
  :param id_column: column index of id in input file
  :param mapping_source_column: column index of source ids in mapping file
  :param mapping_converted_column: column index of converted ids in mapping file
  """
  mappings = parse_mapping(mapping_file=mapping_file,
                           source_column=mapping_source_column,
                           converted_column=mapping_converted_column)

  for line in input_file:
    if line.startswith('#'):
      continue
    columns = line.rstrip('\r\n').split('\t')
    source = columns[id_column]
    mapping = mappings[source] if source in mappings else ''
    if not mapping and source.find('.') >= 0:
      source = source[:source.index('.')]
      mapping = mappings[source] if source in mappings else ''
    new_columns = columns[0:id_column + 1]
    new_columns.append(mapping)
    new_columns.extend(columns[id_column + 1:])
    output_file.write('\t'.join(new_columns))
    output_file.write('\n')


def parse_mapping(mapping_file: TextIO, source_column: int = 0,
    converted_column: int = 1) \
    -> dict[str, str]:
  """
  Parse mapping file.

  :param mapping_file: text delimited filename
  :param source_column: index of source id columns
  :param converted_column: index of converted id columns
  :return: dictionary of source id to converted id
  """
  mappings = {}
  for line in mapping_file:
    if line.startswith('#'):
      continue
    columns = line.rstrip('\r\n').split('\t')
    source = columns[source_column]
    converted = columns[converted_column]
    mappings[source] = converted
  return mappings


if __name__ == '__main__':
  main()
