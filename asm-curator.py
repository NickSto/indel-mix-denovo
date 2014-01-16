#!/usr/bin/env python
import os
import sys
import random
import quicksect
import bioreaders
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'bool':False}
USAGE = "USAGE: %prog [options] lastz-output.lav assembly.fa"
DESCRIPTION = """"""
EPILOG = """"""

#TODO: handle contigs with exactly the same start/end points
#      maybe use IntervalNode.linenum as unique identifier
def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-s', '--str', dest='str',
    default=OPT_DEFAULTS.get('str'), help='default: %default')
  parser.add_option('-i', '--int', dest='int', type='int',
    default=OPT_DEFAULTS.get('int'), help='')
  parser.add_option('-f', '--float', dest='float', type='float',
    default=OPT_DEFAULTS.get('float'), help='')
  parser.add_option('-b', '--bool', dest='bool', action='store_true',
    default=OPT_DEFAULTS.get('bool'),
    help='')

  (options, arguments) = parser.parse_args()

  if len(arguments) != 2:
    parser.print_help()
    fail("\nError: Please provide exactly two positional arguments.")
  else:
    for path in arguments:
      if not os.path.isfile(path):
        fail("\nError: cannot find file "+path)
    (lavpath, fastapath) = arguments

  lav = bioreaders.LavReader(lavpath)
  lav.convert()
  intervals = lav_to_intervals(lav)
  all_overlaps = get_all_overlaps(intervals)
  for interval in sorted(intervals, key=lambda x: x[0]):
    print "overlapping: ",format_interval(interval, intervals)
    for overlap in sorted(all_overlaps[interval], key=lambda x: x[0]):
      print format_interval(overlap, intervals)


def format_interval(interval, intervals):
  output = str(interval)
  output += ' '+str(interval[1]-interval[0]+1)+' bp:\t'
  name = intervals[interval].parent.query['id']
  if name.startswith('NODE_'):
    fields = name.split('_')
    if len(fields) > 2:
      name = '_'.join(fields[:2])
  output += name
  return output


def lav_to_intervals(lav):
  """Convert a LASTZ alignment to a series of intervals along the reference
  sequence.
  Give an LavReader and it will return a dict with the intervals as keys and
  the corresponding LavAlignments as values. Each interval is a tuple of the
  "begin" and "end" values of an LavAlignment."""
  intervals = {}
  for hit in lav:
    for alignment in hit:
      interval = (alignment.subject['begin'], alignment.subject['end'])
      intervals[interval] = alignment
  return intervals


def get_all_overlaps(intervals):
  # build tree
  tree = None
  random.seed(1)
  for interval in intervals:
    if tree is None:
      tree = quicksect.IntervalNode(interval[0], interval[1])
    else:
      tree = tree.insert(interval[0], interval[1])
  # find ones that intersect each interval
  all_overlaps = {}
  for interval in intervals:
    overlaps = []
    # reporter function: is given the matching interval when one is found
    add_overlap = lambda node: overlaps.append((node.start, node.end))
    tree.intersect(interval[0], interval[1], add_overlap)
    # remove the query interval from the results
    overlaps[:] = [overlap for overlap in overlaps if overlap != interval]
    all_overlaps[interval] = overlaps
  return all_overlaps


def intlen(interval):
  return interval[1] - interval[0] + 1


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
