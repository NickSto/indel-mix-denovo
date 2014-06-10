#!/usr/bin/env python
import random
import quicksect
"""Methods useful for manipulating intervals in an LAV file."""
__version__ = '0.7'

SLOP_DEFAULT = 20

def alignments_to_intervals(lav, reference='subject'):
  """Convert a set of LASTZ alignments to a set of intervals along the reference
  sequence.
  Give an LavReader and it will return a dict with the intervals as keys and
  the corresponding LavAlignments as values. Each interval is a tuple of the
  "begin" and "end" values of an LavAlignment.
  The reference sequence is assumed to be the "subject", but change "reference"
  to "query" to use the other sequence.
  N.B.: Alignments with identical begin/end coordinates are not distinguished
  (only the last one encountered will be kept)."""
  assert reference in ('subject', 'query'), (
    '"reference" must be "subject" or "query".')
  intervals = {}
  for hit in lav:
    for alignment in hit:
      seq = getattr(alignment, reference)
      interval = (seq['begin'], seq['end'])
      intervals[interval] = alignment
  return intervals


def blocks_to_intervals(lav, reference='subject'):
  """Convert a set of LASTZ blocks to a set of intervals along the reference
  sequence.
  Give an LavReader and it will return a dict with the intervals as keys and
  the corresponding LavBlocks as values. Each interval is a tuple of the
  "begin" and "end" values of an LavBlocks.
  The reference sequence is assumed to be the "subject", but change "reference"
  to "query" to use the other sequence.
  N.B.: Blocks with identical begin/end coordinates are not distinguished
  (only the last one encountered will be kept)."""
  assert reference in ('subject', 'query'), (
    '"reference" must be "subject" or "query".')
  intervals = {}
  for hit in lav:
    for alignment in hit:
      for block in alignment:
        seq = getattr(block, ref)
        interval = (seq['begin'], seq['end'])
        intervals[interval] = block
  return intervals


def map_subject_to_query(lav, level='alignments'):
  """Convert LASTZ alignments to intervals along both the subject and query
  sequences, and return a dict mapping between the two.
  Each key of the dict will be the interval of one alignment along the subject
  sequence, and its value will be the alignment's interval along the query."""
  if level != 'alignments':
    raise NotImplementedError
  intervals = {}
  for hit in lav:
    for alignment in hit:
      subject_interval = (alignment.subject['begin'], alignment.subject['end'])
      query_interval = (alignment.query['begin'], alignment.query['end'])
      intervals[subject_interval] = query_interval
  return intervals


def conversion_coefficients(block):
  """Return the coefficients needed to convert a query coordinate in "block"
  into a subject coordinate.
  Returns a tuple of (strand, offset), where
  subject_coord = query_coord * strand + offset
  strand is -1 if the query is reverse complemented, 1 otherwise
  """
  if block.parent.parent.query['revcomp']:
    strand = -1
  else:
    strand = 1
  offset = block.subject['begin'] - strand * block.query['begin']
  return (strand, offset)


def blocks_to_conv_table(lav, contigs=None):
  """Build a coordinate conversion table from an LAV file.
  The table lists the values needed to convert sites in each gap-free block
  from query coordinates to subject coordinates.
  "lav" is an LavReader object.
  "contigs" contains the names of the valid query sequences to add to the table.
    If "contigs" is given, any hit whose query name is not in "contigs" will be
    left out of the table.
  The return value is a list of tuples, one per LavBlock. The 5 elements of
  each tuple are:
  0 (chrom):  the chromosome the block is in (the query name)
  1 (begin):  the block's start coordinate (in the query)
  2 (end):    the block's end coordinate (in the query)
  3 (ref):    the name of the corresponding reference chromosome
  4 (strand): the strand value for coordinate conversion
  5 (offset): the offset for coordinate conversion
  The last two values are the output of conversion_coefficients().
  """
  intervals = blocks_to_intervals(lav)
  table = []
  for hit in lav:
    chrom = hit.query['name'].split()[0]
    ref = hit.subject['name'].split()[0]
    if contigs is not None and chrom not in contigs:
      continue
    for alignment in hit:
      for block in alignment:
        if block.query['begin'] < block.query['end']:
          begin = block.query['begin']
          end = block.query['end']
        else:
          begin = block.query['end']
          end = block.query['begin']
        coef = conversion_coefficients(block)
        table.append([chrom, begin, end, ref, coef[0], coef[1]])
  return table


def get_all_overlaps(intervals, sort=False):
  """Compare each interval to the rest, finding ones with any overlap.
  Returns a dict mapping each interval to a list of intervals which overlap.
  The original interval is always excluded from the list.
  When "sort" is True, each list of overlapping intervals will be sorted
  according to starting coordinate."""
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
    if sort:
      overlaps.sort(key=lambda interval: interval[0])
    all_overlaps[interval] = overlaps
  return all_overlaps


def get_overlaps(tree, interval):
  """Return a list of the intervals overlapping with a query interval.
  "interval" is the query interval, a tuple of (start, stop) coordinates.
  "tree" is a quicksect.IntervalNode of set of intervals which might overlap.
  If "tree" is None, this will return an empty list."""
  overlaps = []
  if tree is None:
    return overlaps
  # reporter function: is given the matching interval when one is found
  add_overlap = lambda node: overlaps.append((node.start, node.end))
  tree.intersect(interval[0], interval[1], add_overlap)
  return overlaps


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