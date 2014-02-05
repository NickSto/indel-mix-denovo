#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import lavreader
import lavintervals

OPT_DEFAULTS = {'float':1.0}
USAGE = "USAGE: %(prog)s [options] align.lav sites.tsv"
DESCRIPTION = """Default input format: the tsv output of inspect-reads.py."""
EPILOG = """"""

SLOP = 20

def main():

  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('lavpath', metavar='align.lav',
    help="""""")
  parser.add_argument('sitespath', metavar='sites.tsv',
    help="""""")
  parser.add_argument('-v', '--vcf', action='store_true',
    help="""Input data (second file) is a VCF.""")
  parser.add_argument('-f', '--float', type=float,
    help="""default: %(default)s""")

  args = parser.parse_args()

  lav = lavreader.LavReader(args.lavpath)
  contigs = get_kept_contigs(lav)
  table = blocks_to_conv_table(lav, contigs=contigs)


  # for each site:
  # just compare to each interval to check if it's contained
  # use the longest one


def get_kept_contigs(lav):
  contigs = set()
  #TODO: replace with reading retained contigs from FASTA
  intervals = lavintervals.alignments_to_intervals(lav)
  all_overlaps = lavintervals.get_all_overlaps(intervals)
  all_overlaps = lavintervals.discard_redundant(all_overlaps, slop=SLOP)
  # remove discarded intervals
  for interval in intervals.keys():
    if interval in all_overlaps:
      hit = intervals[interval].parent
      contigs.add(hit.query['name'])
  return contigs


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
