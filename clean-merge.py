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
    help='Input LAV file aligning a single query sequence to a single '
      'reference (subject) sequence. E.g. from a command like "$ lastz chrM.fa '
      'assembly.fa > align.lav".')
  parser.add_argument('fasta', metavar='assembly.fa',
    help='The assembly FASTA file (the query in the LAV alignment).')
  parser.add_argument('-b', '--bool', action='store_true',
    help='(stores False by default if not given)')

  args = parser.parse_args()

  # open LAV and FASTA files
  fasta = fastareader.FastaLineGenerator(args.fasta)
  lav = lavreader.LavReader(args.lav)
  if len(lav) > 1:
    fail('Error: Found more than 1 subject and/or query sequence in alignment.')

  # Compute info that will be needed on the intervals:
  # convert alignments to a set of intervals and map back to lav objects
  interval_to_aln = lavintervals.alignments_to_intervals(lav)
  # extract the intervals into a simple list
  intervals = interval_to_aln.keys()
  # compute which of the other intervals overlap each one
  all_overlaps = lavintervals.get_all_overlaps(intervals)
  # build a table to convert subject coordinates to query coordinates
  conversion_table = lavintervals.blocks_to_conv_table(lav, dicts=True,
    query_to_subject=False)

  # Construct a final, non-redundant sequence out of the original, by walking
  # along the reference, adding sequence from the assembly.
  # Algorithm:
  # Sort list of intervals by starting coordinate
  # While intervals left:
  #   Pop the first interval off the sorted list
  #   If the interval is unique (does not overlap any other interval):
  #     Take the corresponding sequence from the assembly & add it to the output
  #   Else (it overlaps):
  #     Break the interval in two: a unique region, followed by the rest
  #     Add the unique region to the list of intervals
  #     
  # 
  # Walk along the reference, adding contigs when they uniquely cover a region,
  # and if there are multiple contigs to choose from, decide using LASTZ scores
  # of alignments to the original assembly.
  # Uses "intervals" as an ordered queue.
  final_sequence = ''
  while len(intervals) > 1:
    # sort by start coord
    #TODO: add new intervals intelligently, to preserve order & avoid sorting
    intervals.sort(key=lambda interval: interval[0])
    interval = intervals[0]
    next_interval = intervals[1]
    intervals.remove(interval)
    # Is this interval unique (no overlap between this interval and the next)?
    if interval[1] < next_interval[0]:
      # Add query sequence of the interval
      aln = interval_to_aln[interval]
      final_sequence += fasta.extract(aln.query['begin'], aln.query['end'])
      # Is there a gap between this interval and the next?
      if interval[1] - next_interval[0] > 1:
        pass #TODO: add N's? create a new sequence?
    # or do the intervals overlap?
    else:
      unique = (interval[0], next_interval[0]-1)
      rest = (next_interval[0], interval[1])
      overlap = (next_interval[0], min(interval[1], next_interval[1]))
      # Add the unique sequence before the overlap (if any) to the list
      if unique[0] <= unique[1]:
        unique_on_query = convert_interval(conversion_table, unique)
        final_sequence += fasta.extract(*unique_on_query)
      # Determine which sequence to use in the overlap section
      if choose_contig(interval, next_interval, fasta, interval_to_aln):
        pass
      else:
        pass
    # break up into new intervals if needed, add entries in interval_to_aln
  # add final interval


  # New idea: extract FASTA for each redundant hit, LASTZ align to original
  # assembly, and choose the one with the highest-scoring hit in each case.


def get_query_seq(subinterval, interval, interval_to_aln, fasta):
  """"""

def convert_interval(table, interval):
  begin = lavintervals.convert(table, interval[0])
  end = lavintervals.convert(table, interval[1])
  return (begin, end)

def interval_to_fasta(intervals, fasta, fastapath):
  """Extract the given interval from a fasta file and write the sequence into
  a new file."""

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
