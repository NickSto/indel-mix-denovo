#!/usr/bin/env python
from __future__ import division
import os
import sys
import random
import argparse
import lavreader
import fastareader
import lavintervals

OPT_DEFAULTS = {'str':'string', 'int':0}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
EPILOG = """"""

def main():

  parser = argparse.ArgumentParser(
    description=DESCRIPTION, usage=USAGE, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('lav', metavar='align.lav',
    help='Input LAV file. Must be an alignment of a single query sequence to a '
      'single subject sequence.')
  parser.add_argument('fasta', metavar='assembly.fa',
    help='The assembly FASTA file (the query in the LAV alignment).')
  parser.add_argument('-b', '--bool', action='store_true',
    help='(stores False by default if not given)')

  args = parser.parse_args()

  # open, parse LAV and FASTA files
  fasta = fastareader.FastaLineGenerator(args.fasta)
  lav = lavreader.LavReader(args.lav)
  if len(lav) > 1:
    fail('Error: Found more than 1 subject and/or query sequence in alignment.')
  # convert alignments to a set of intervals
  intervals = lavintervals.alignments_to_intervals(lav)
  # get list of other alignments that overlap each alignment
  all_overlaps = lavintervals.get_all_overlaps(intervals)

  # New idea: extract FASTA for each redundant hit, LASTZ align to original
  # assembly, and choose the one with the highest-scoring hit in each case.

  # concatenate alignments separated only by Ns
  subj_to_query = lavintervals.map_subject_to_query(lav)
  for (interval, overlaps) in all_overlaps.items():
    all_overlaps[interval] = join_intervals(overlaps, fasta)


def interval_to_fasta(intervals, fasta, fastapath):
  """Extract the given interval from a fasta file and write the sequence into
  a new file."""

def get_seq(interval, fasta):
  """Give a (begin, end) interval and a FastaLineGenerator."""
  return fasta.extract(*interval)

def get_tmp_path(base):
  attempts = 1
  candidate = base + '.tmp'
  while os.path.exists(candidate):
    candidate = "{}.{}.{}".format(base, random.randint(0,1000), 'tmp')
    attempts+=1
    if attempts > 20:
      fail('Error: cannot find an unoccupied temporary directory name.')
  return candidate


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
