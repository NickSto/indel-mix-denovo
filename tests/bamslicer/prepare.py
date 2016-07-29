#!/usr/bin/env python2
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import argparse
import subprocess

ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Prepare test input files for bamslicer.py from a human-readable file of reads and
variants."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('infile')
  parser.add_argument('ref',
    help='Output reference file (fasta).')
  parser.add_argument('reads',
    help='Output reads file (fastq).')
  parser.add_argument('expected',
    help='Expected output file (tsv) which bamslicer will write to $TEST_OUTPUT.')
  parser.add_argument('-a', '--align', action='store_true',
    help='Also align the reads to the reference using bowtie2 (must be on your $PATH).')

  args = parser.parse_args(argv[1:])

  # Parse the input file.
  variants = []
  read_locations = []
  with open(args.infile) as lines, \
       open(args.ref, 'w') as ref_file, \
       open(args.reads, 'w') as reads_file:
    for line in lines:
      # comment
      if line.strip().startswith('#'):
        continue
      # reference
      if line.startswith('r'):
        ref_file.write('>ref\n')
        ref_file.write(line[1:])
      # reads
      if line.startswith('R'):
        matches = re.finditer(r'(\d+):([ACGT]+)', line)
        for match in matches:
          name = match.group(1)
          seq  = match.group(2)
          reads_file.write('@{}\n'.format(name))
          reads_file.write(seq+'\n')
          reads_file.write('+\n')
          reads_file.write('H' * len(seq) + '\n')
          read_locations.append((name, match.start(2), match.end(2)-1))
      # variants
      if line.startswith('v'):
        coord = -1
        for char in line:
          coord += 1
          if coord == 0:
            continue
          if char in ('A', 'G', 'T', 'C'):
            variants.append((coord, char))

  # Print the expected output from bamslicer.py's debugging.
  with open(args.expected, 'w') as expected:
    for coord, base in variants:
      names = []
      for name, start, end in read_locations:
        if start <= coord <= end:
          names.append(name)
      if names:
        expected.write('{}\t{}\t{}\t{}\n'.format(coord, 'I', base, ','.join(names)))

  # Print the variants, in notation understood by inspect-reads.py.
  print(','.join(['ref:{}-I:{}'.format(coord, char) for coord, char in variants]))
  sys.stderr.write('\n')

  if args.align:
    align(args.ref, args.reads)


def align(ref, reads):
  outbase = os.path.splitext(reads)[0]
  if outbase == reads:
    outbase = reads + '.out'
  with open(os.devnull, 'w') as devnull:
    subprocess.call(['bowtie2-build', ref, ref], stdout=devnull)
  subprocess.call(['bowtie2', '-x', ref, '-U', reads, '-S', outbase+'.sam'])
  with open(outbase+'.tmp.bam', 'w') as bamtmp, open(os.devnull, 'w') as devnull:
    subprocess.call(['samtools', 'view', '-Sb', outbase+'.sam'], stdout=bamtmp, stderr=devnull)
  subprocess.call(['samtools', 'sort', outbase+'.tmp.bam', outbase])
  subprocess.call(['samtools', 'index', outbase+'.bam'])
  sys.stderr.write('The aligned bam file is "'+outbase+'.bam"\n')


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
