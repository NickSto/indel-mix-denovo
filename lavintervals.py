#!/usr/bin/env python
"""Methods useful for manipulating intervals in an LAV file."""
import random
import quicksect

SLOP_DEFAULT = 20

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
    overlaps = [overlap for overlap in overlaps if overlap != interval]
    all_overlaps[interval] = overlaps
  return all_overlaps


def discard_redundant(all_overlaps, slop=SLOP_DEFAULT):
  """Remove redundant intervals from all_overlaps.
  The discarded intervals will be removed from both the
  set of keys and the lists of overlapping intervals."""
  redundant = set()
  # gather redundant intervals
  for interval in all_overlaps:
    if is_redundant(interval, all_overlaps[interval], slop=slop):
      redundant.add(interval)
  # go through all_overlaps again, discarding redundant intervals
  for interval in all_overlaps.keys():
    if interval in redundant:
      del(all_overlaps[interval])
    else:
      all_overlaps[interval] = [overlap for overlap in all_overlaps[interval]
        if overlap not in redundant]
  return all_overlaps


def is_redundant(interval, overlaps, slop=SLOP_DEFAULT):
  """Return True if the interval is redundant."""
  for overlap in overlaps:
    # if interval is wholly contained within overlap (+ slop)
    if (overlap[0] - slop <= interval[0] and interval[1] <= overlap[1] + slop
        and length(interval) < length(overlap)):
      return True
  return False


def length(interval):
  return interval[1] - interval[0] + 1


def merge(intervals, sort=False):
  """Merge a list of intervals into a new list of non-overlapping intervals.
  Optionally sort the output list by starting coordinate. Merging is 1-based
  (will merge if end == start). Totally naive implementation: O(n^2)."""
  merged = []
  # go through input list, adding 1 interval at a time to a growing merged list
  for interval in intervals:
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
  if sort:
    merged.sort(key=lambda x: x[0])
  return merged


def subtract(merged, interval):
  """Subtract intervals in merged from interval. Merged must be a list of
  "merged" intervals: non-overlapping, in order from lowest to highest
  coordinate. Returns a list of intervals representing the regions of interval
  not overlapping any intervals in merged."""
  begin = interval[0]
  end = interval[1]
  results = []
  for i in range(len(merged)):
    # add region before first overlap (if any)
    if i == 0 and begin < merged[i][0]:
      results.append((begin, merged[i][0]-1))
    if i > 0:
      result = (merged[i-1][1]+1, merged[i][0]-1)
      if result[0] <= result[1]:
        results.append(result)
    # add region after last overlap (if any)
    if i == len(merged) - 1 and end > merged[i][1]:
      results.append((merged[i][1]+1, end))
  if not merged:
    results.append((begin, end))
  return results