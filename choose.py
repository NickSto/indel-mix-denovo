#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import argparse

ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Quick script to choose assemblies from families based on directory name and number
of contigs."""
EPILOG = """Assumes directory names like M117-bl, M117-ch, M117C1-bl, M117C1-ch"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('asm_dir_root',
    help='The directory which contains the SPAdes output directories.')

  args = parser.parse_args(argv[1:])

  asm_dirs = get_asm_dirs(args.asm_dir_root)

  families = group_dirs(asm_dirs)

  for family in families:
    min_contigs = 99999999
    min_asm = None
    for asm_dir in families[family]:
      asm_path = os.path.join(asm_dir, 'contigs.fasta')
      contigs = count_contigs(asm_path)
      if contigs > 0 and contigs < min_contigs:
        min_contigs = contigs
        min_asm = os.path.basename(asm_dir)
    print(family, min_asm, min_contigs, sep='\t')


def get_asm_dirs(asm_dir_root):
  asm_dirs = []
  for name in os.listdir(asm_dir_root):
    path = os.path.join(asm_dir_root, name)
    if os.path.isdir(path):
      asm_dirs.append(path)
  return sorted(asm_dirs)


def group_dirs(asm_dirs):
  families = {}
  for asm_dir in asm_dirs:
    dirname = os.path.basename(asm_dir)
    match = re.search(r'^([A-Za-z]+\d+)(C[0-9]+)?', dirname)
    if match:
      family = match.group(1)
      samples = families.get(family, [])
      samples.append(asm_dir)
      families[family] = samples
  return families


def count_contigs(asm_path):
  if not os.path.isfile(asm_path):
    return 0
  contigs = 0
  with open(asm_path) as asm:
    for line in asm:
      if line.startswith('>'):
        contigs += 1
  return contigs


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
