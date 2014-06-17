#!/usr/bin/env python
from __future__ import division
import os
import sys
import random
import shutil
import argparse
import subprocess
import lavintervals
import fastareader
import lavreader

OPT_DEFAULTS = {'str':'string', 'int':0}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
EPILOG = """"""
FASTA_WIDTH = 70
TMP_DIR_BASE = 'cleaning'

def main():

  parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('lav', metavar='align.lav',
    help='Input LAV file aligning a single query sequence to a single '
      'reference (subject) sequence. E.g. from a command like "$ lastz chrM.fa '
      'assembly.fa > align.lav".')
  parser.add_argument('asm_merge', metavar='assembly.fa',
    help='The FASTA file of the assembly to be cleaned (the query in the LAV '
      'alignment).')
  parser.add_argument('asm_raw', metavar='raw-assembly.fa',
    help='The FASTA file of the original assembly, before curating, merging, '
      'etc.')

  args = parser.parse_args()

  # open LAV and FASTA files, make temporary directory
  asm_merge = fastareader.FastaLineGenerator(args.asm_merge)
  lav = lavreader.LavReader(args.lav)
  if len(lav) > 1:
    fail('Error: Found more than 1 subject and/or query sequence in alignment.')
  tmpdir = get_tmp_path(TMP_DIR_BASE)
  try:
    os.makedirs(tmpdir)
  except OSError:
    fail('Error: temporary directory "'+tmpdir+'" exists.')


  # Compute info that will be needed on the intervals:
  # convert alignments to a set of intervals and map back to lav objects
  interval_to_aln = lavintervals.alignments_to_intervals(lav)
  # extract the intervals into a simple list
  intervals = interval_to_aln.keys()
  # build a table to convert subject coordinates to query coordinates
  conversion_table = lavintervals.hits_to_conv_table(lav, query_to_subject=False)

  # Construct a final, non-redundant sequence out of the original, by walking
  # along the reference, adding sequence from the assembly.
  # Algorithm:
  # While intervals left:
  #   Sort list of intervals by starting coordinate
  #   Pop the first interval off the list
  #   If the interval is unique (does not overlap any other interval):
  #     Take the corresponding sequence from the assembly & add it to the output
  #   Else (it overlaps):
  #     Break the interval in two: a unique region, followed by the rest
  #     Add the unique region (if any) to the list of intervals
  #     Decide between the two overlapping intervals:
  #       For each, align the overlapping region to the original assembly
  #       Use the one with the highest LASTZ score
  #     If the larger interval wins:
  #       Pop the smaller interval from the list
  #     Else (smaller interval wins):
  #       Split the larger interval and pop it from the list
  #       Add the non-overlapping segment to the list
  final_sequence = ''
  while len(intervals) > 0:
    # Pop the next interval, peek the one after it (if it exists)
    #TODO 3: add new intervals intelligently, to preserve order & avoid sorting
    intervals.sort(key=lambda interval: interval[0]) # sort by start coord
    interval = intervals[0]
    if len(intervals) > 1:
      next_interval = intervals[1]
    else:
      next_interval = None
    del(intervals[0])
    assert next_interval is None or interval[0] <= next_interval[0]
    # Is this interval unique (no overlap between this interval and the next)?
    if next_interval is None or interval[1] < next_interval[0]:
      # Add query sequence of the interval to the output sequence
      interval_on_query = convert_with_alignment(interval,
                                                 interval_to_aln[interval])
      final_sequence += asm_merge.extract(*interval_on_query)
      # Is there a gap between this interval and the next?
      if next_interval is not None and interval[1] - next_interval[0] > 1:
        #TODO 2: add N's (goal is to keep reference coordinates)
        #        not necessary for R33S10 (those gaps are in redundant seq)
        raise NotImplementedError('There must be no gaps between contigs.')
    # or do the intervals overlap?
    else:
      # The overlapping region
      overlap = (next_interval[0], min(interval[1], next_interval[1]))
      # Portion of interval before the overlapping region
      unique = (interval[0], next_interval[0]-1)
      # All parts of interval after the unique region
      rest = (next_interval[0], interval[1])
      # Add the unique sequence before the overlap (if any) to the list
      if length(unique) > 0:
        intervals.append(unique)
        interval_to_aln[unique] = interval_to_aln[interval]
      # Determine which sequence to use in the overlap section
      alignment1 = interval_to_aln[interval]
      alignment2 = interval_to_aln[next_interval]
      if choose_sequence(alignment1, alignment2, overlap, asm_merge,
                         args.asm_raw, tmpdir):
        (winner, loser) = (rest, next_interval)
      else:
        (winner, loser) = (next_interval, rest)
      assert winner[0] == loser[0]
      # Make sure the list contains the winning interval and not the losing one.
      if loser in intervals:
        intervals.remove(loser)
      if winner not in intervals:
        intervals.append(winner)
        if winner == next_interval:
          interval_to_aln[winner] = interval_to_aln[next_interval]
        else:
          interval_to_aln[winner] = interval_to_aln[interval]
      # If the losing interval is the larger one, split it and keep the part
      # after the overlap.
      if length(loser) > length(winner):
        loser_rest = (winner[1]+1, loser[1])
        intervals.append(loser_rest)
        if loser == next_interval:
          interval_to_aln[loser_rest] = interval_to_aln[next_interval]
        else:
          interval_to_aln[loser_rest] = interval_to_aln[interval]
  print fasta_format(final_sequence, 'Cleaned', FASTA_WIDTH)

  shutil.rmtree(tmpdir)


def length(interval):
  """Get length of interval.
  1-based: length((10, 10)) == 1"""
  return interval[1] - interval[0] + 1


def choose_sequence(alignment1, alignment2, overlap, asm_merge, asm_raw_file,
                    tmpdir):
  """Returns True if the first sequence is best, False otherwise."""
  top_scores = []
  for (alignment, name) in zip((alignment1, alignment2), ('seq1', 'seq2')):
    fastapath = os.path.join(tmpdir, name+'.fa')
    overlap_on_query = convert_with_alignment(overlap, alignment)
    sys.stderr.write("{} -> {}\n".format(overlap, overlap_on_query))
    sequence = asm_merge.extract(*overlap_on_query)
    with open(fastapath, 'w') as fastafile:
      fastafile.write(fasta_format(sequence, name))
    # Perform LASTZ alignment
    lavpath = os.path.join(tmpdir, name+'.lav')
    with open(lavpath, 'w') as lavfile:
      # sys.stderr.write("$ "+" ".join(['lastz', fastapath, asm_raw_file])+"\n")
      subprocess.call(['lastz', fastapath, asm_raw_file], stdout=lavfile)
    lav = lavreader.LavReader(lavpath)
    # Get the top-scoring alignment
    top_score = 0
    for hit in lav:
      for alignment in hit:
        if alignment.score > top_score:
          top_score = alignment.score
    top_scores.append(top_score)
  if top_scores[0] > top_scores[1]:
    return True
  else:
    return False


def fasta_format(sequence, name, width=FASTA_WIDTH):
  """Turn a sequence and name into a FASTA-formatted string."""
  output = '>'+name+'\n'
  for start in range(0, len(sequence), width):
    end = start + width
    if end > len(sequence):
      end = len(sequence)
    output += sequence[start:end]+'\n'
  return output


def convert_with_alignment(interval, alignment):
  """One-off interval conversion, using a pre-chosen alignment.
  Converts the interval's start/end coordinates from subject to query
  coordinates using the given alignment.
  Assumes the interval is actually contained in the alignment."""
  #TODO 1: What happens when a coordinate is in a small gap between blocks?
  #        Make sure the resulting sequence doesn't include 1bp duplications.
  table = lavintervals.alignments_to_conv_table([alignment],
    query_to_subject=False)
  try:
    begin = lavintervals.convert(table, interval[0], fail='throw')[0]
    end = lavintervals.convert(table, interval[1], fail='throw')[0]
  except Exception as e:
    if len(e.args) > 1 and e.args[1] == 'fail':
      raise AssertionError('Interval must be contained in alignment.')
    else:
      raise
  if begin > end:
    (begin, end) = (end, begin)
  return (begin, end)


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
