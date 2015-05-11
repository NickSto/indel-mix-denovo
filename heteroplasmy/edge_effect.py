#!/usr/bin/env python
# Input: File with 4 columns: sample, position, major, minor allele
# also need a BAM file named by the sample name from the input file
# - uses read groups from the BAM file to get the sample name
from __future__ import division
import os
import sys
import pysam
import argparse
import distutils.version

OPT_DEFAULTS = {'chrom':'chrM', 'min_qual':30, 'edge_size':25}

USAGE = "$ %(prog)s heteroplasmies.tsv reads.bam [outfile.tsv]"
DESCRIPTION = """Determine whether a variant is mostly supported by reads where
the variant appears in edge of the read."""


def main():

  # Broken by PySAM 0.8.1.
  # Fails with error:
  # AttributeError: 'pysam.calignmentfile.PileupRead' object has no attribute 'qpos'
  # on line:
  #             read.alignment.seq[read.qpos] == minor and
  # Looks they probably renamed qpos to something like query_alignment_position, though they don't
  # specifically mention it in the release notes:
  # http://pysam.readthedocs.org/en/latest/release.html#release-0-8-1
  version = pysam.__version__
  if distutils.version.StrictVersion(version) >= distutils.version.StrictVersion('0.8.1'):
    raise Exception('Not compatible with PySAM versions > 0.8.0. Current version: '+version)

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('varfile', metavar='heteroplasmies.tsv',
    help='Input variants. If -v is not given, should have 4 columns: sample, '
      'position, major, minor allele.')
  parser.add_argument('bamfile', metavar='reads.bam',
    help='Alignment the variants were called in. Must use read groups with the '
      'same sample names as in the variants tsv file.')
  parser.add_argument('outfile', metavar='outfile.tsv', nargs='?',
    help='The output file. If not given, will print to stdout. Output is '
      'tab-separated, with 5 columns: 1. sample name. 2. coordinate. 3. number '
      'of reads where the minor allele appears in a read edge. 4. total number '
      'of minor allele reads. 5. ratio of column 3/4')
  parser.add_argument('-v', '--va-input', action='store_true',
    help='Give if the variants tsv file is the raw output of the Variant '
      'Annotator. Must be unstranded output with either 12 or 13 columns '
      '(strand bias column 13 is optional).')
  parser.add_argument('-a', '--append', action='store_true',
    help='Append the edge statistics to the input file. Instead of printing '
      'just the 5 columns mentioned above, print each input line, with 3 '
      'additional columns (the last 3 of the normal output).')
  parser.add_argument('-c', '--chrom',
    help='The chromosome to restrict analysis to. Default: %(default)s.')
  parser.add_argument('-q', '--min-qual', type=int,
    help='Minimum base quality. If the variant is present in a read, but its '
      'base call lower than this quality, it won\'t be counted (as either in '
      'the edge or not). Default: %(default)s.')
  parser.add_argument('-e', '--edge-size', type=int,
    help='Size, in bases, of the read edge region. Default: %(default)s.')

  args = parser.parse_args()

  if args.outfile:
    outfile = open(args.outfile, 'w')
  else:
    outfile = sys.stdout

  if args.va_input:
    (sample_col, pos_col, major_col, minor_col) = (0, 2, 9, 10)
  else:
    (sample_col, pos_col, major_col, minor_col) = (0, 1, 2, 3)

  with open(args.varfile) as data:
    for line in data:
      # skip empty lines
      if not line:
        continue
      fields = line.strip().split('\t')
      sample = fields[sample_col]
      pos = int(fields[pos_col])
      major = fields[major_col]
      minor = fields[minor_col]
      if sample == '__NONE__':
        fail('Error: sample name is __NONE__')
      edge_effect = get_edge_effect(args, sample, pos, major, minor)
      if args.append:
        outfile.write(line.rstrip('\r\n')+'\t'+'\t'.join(edge_effect)+'\n')
      else:
        outfile.write('\t'.join([sample, str(pos)] + edge_effect)+'\n')

  if outfile is not sys.stdout:
    outfile.close()


def get_edge_effect(args, sample, pos, major, minor):
  tot_minor_reads = 0
  minor_in_edge = 0
  sam = pysam.Samfile(args.bamfile, 'rb')
  pileup = sam.pileup(args.chrom, pos-1, pos, stepper='all', max_depth=10000000,
                      mask=False, truncate=True)
  for pileupcolumn in pileup:
    if pileupcolumn.pos == pos-1:
      for read in pileupcolumn.pileups:
        # Get the read group
        for (tag, value) in read.alignment.tags:
          if tag == 'RG':
            read_group = value
            break
        # Select reads from the right read group, with the minimum quality,
        # and with the minor allele
        if (read_group == sample and
            read.alignment.seq[read.qpos] == minor and
            ord(read.alignment.qual[read.qpos])-33 >= args.min_qual):
          tot_minor_reads+=1
          rlen = read.alignment.rlen
          edges = range(0, args.edge_size) + range(rlen-args.edge_size, rlen+1)
          if read.qpos in edges:
            minor_in_edge+=1
  try:
    ratio = minor_in_edge/tot_minor_reads
    ratio_str = '%.3f' % (ratio)
  except ZeroDivisionError:
    ratio_str = 'na'
  return [str(minor_in_edge), str(tot_minor_reads), ratio_str]


def fail(message):
  sys.stderr.write(message+'\n')
  sys.exit(1)


if __name__ == '__main__':
  main()
