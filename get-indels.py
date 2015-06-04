#!/usr/bin/env python
from __future__ import division
import sys
import argparse
import collections
import pyBamParser.bam

OPT_DEFAULTS = {'threshold':0}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Tally all the indels in a BAM alignment and print their read counts. Prints output
to stdout: A list of indels and their counts, one indel per line. Each line is 4 tab-delimited
columns: coordinate, indel type (D or I), insertion or deletion length, and count."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('bam', metavar='align.bam',
    help='Your read alignment.')
  parser.add_argument('-t', '--threshold', type=int,
    help='The minimum number of reads that need to support an indel for it to be reported. '
         'Default: %(default)s')
  # parser.add_argument('-r', '--region',
  #   help='Restrict to this region. Give in UCSC format, e.g. "chr1:500-600".')
  parser.add_argument('-g', '--read-group',
    help='Only look at this read group.')
  # parser.add_argument('-G', '--read-groups', action='store_true',
  #   help='Separate counts by read group.')

  args = parser.parse_args(argv[1:])

  indels = count_indels(args.bam, read_group=args.read_group)

  print_indels(indels, args.threshold)


def count_indels(bamfile, read_group=None):
  """Examine all indels in all reads, and tally identical indels.
  Returns a dict mapping indels to counts. Each indel is identified with a tuple of
  (position, type, length), where position is the indel coordinate (int), type is 'D' for deletion
  or 'I' for insertion, and length is the deletion or insertion length (int)."""
  bam_reader = pyBamParser.bam.Reader(bamfile)
  indels = collections.defaultdict(int)
  for read in bam_reader:
    sample = read.get_read_group()
    # Skip if the user specified a read group and this isn't it.
    if read_group and sample != read_group:
      continue
    # read_pos = read.get_position()
    # read_end = read.get_end_position()
    (read_ins, read_del) = read.get_indels()
    for (del_pos, del_len) in read_del:
      indel = (del_pos, 'D', del_len)
      indels[indel] += 1
    for (ins_pos, ins_len) in read_ins:
      indel = (ins_pos, 'I', ins_len)
      indels[indel] += 1
  return indels


def print_indels(indels, threshold):
  """Print all indel counts above the threshold in a tab-delimited format."""
  for indel in sorted(indels):
    count = indels[indel]
    if count < threshold:
      continue
    (pos, type_, length) = indel
    print "\t".join(map(str, (pos, type_, length, count)))


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
