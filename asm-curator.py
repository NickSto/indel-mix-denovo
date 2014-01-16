#!/usr/bin/env python
import os
import sys
import random
import quicksect
import bioreaders
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'test_output':False}
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
  parser.add_option('-T', '--test-output', dest='test_output',
    action='store_true', default=OPT_DEFAULTS.get('test_output'),
    help='Print legacy test output.')

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
  # all_overlaps = discard_redundant(all_overlaps)
  if options.test_output:
    test_output(all_overlaps, intervals)
    sys.exit()


def test_output(all_overlaps, intervals):
  for interval in sorted(intervals, key=lambda x: x[0]):
    if interval not in all_overlaps:
      continue
    print "overlapping: ",format_interval(interval, intervals)
    for overlap in sorted(all_overlaps[interval], key=lambda x: x[0]):
      print format_interval(overlap, intervals)
    merged = get_unique(all_overlaps[interval])
    print "  merged:"
    for minterval in sorted(merged, key=lambda x: x[0]):
      print "    "+format_interval(minterval)


def format_interval(interval, intervals=None):
  output = str(interval)
  output += ' '+str(interval[1]-interval[0]+1)+' bp'
  if intervals:
    name = intervals[interval].parent.query['id']
    if name.startswith('NODE_'):
      fields = name.split('_')
      if len(fields) > 2:
        name = '_'.join(fields[:2])
    output += ":\t"+name
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
  """Compare each interval to the rest, finding ones with any overlap.
  Returns a dict mapping each interval to a list of intervals which overlap.
  The original interval is always excluded from the list."""
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


def discard_redundant(all_overlaps):
  """Remove intervals wholly contained within other intervals.
  The discarded intervals will be removed from both the set of keys and the
  lists of overlapping intervals."""
  redundant = set()
  # gather redundant intervals
  for interval in all_overlaps:
    if is_redundant(interval, all_overlaps[interval]):
      redundant.add(interval)
    # redundant = redundant.union(get_redundant(interval, all_overlaps))
  # discard redundant intervals
  for interval in all_overlaps.keys():
    if interval in redundant:
      del(all_overlaps[interval])
    else:
      all_overlaps[interval] = filter(lambda x: x not in redundant, all_overlaps[interval])
  return all_overlaps


def is_redundant(interval, overlaps, threshold=0):
  """Return True if the interval is redundant.
  The interval is redundant unless more than threshold % of its length uniquely
  covers the reference."""
  for overlap in overlaps:
    # if interval (the key) is wholly contained within overlap
    if overlap[0] <= interval[0] and interval[1] <= overlap[1]:
      return True
  unique = get_unique(overlaps)


def get_unique(overlaps):
  """Find all regions in overlaps that do not overlap with any other interval.
  Implementation is totally naive, and O(n^2)."""
  merged = []
  for interval in overlaps:
    merged = add_and_merge(interval, merged)
  return merged


def add_and_merge(interval, merged):
  """Add an interval to a list of intervals, merging in the case of overlaps.
  intervals is the growing list of merged intervals. Merging is 1-based
  (will merge if end == start)."""
  #TODO: do this with a sort, then reduce?
  #TODO: then replace with an actually efficient merge, e.g. from bx-python
  i = 0
  while i < len(merged):
    target = merged[i]
    # if any overlap
    if interval[0] <= target[1] and interval[1] >= target[0]:
      # replace the target with a combination of it and the input interval
      del(merged[i])
      # merge the two intervals by taking the widest possible region
      interval = (min(interval[0], target[0]), max(interval[1], target[1]))
    else:
      i+=1
  merged.append(interval)
  return merged


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
