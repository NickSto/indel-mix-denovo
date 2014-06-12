#!/usr/bin/env python
from __future__ import division
import random
import quicksect
"""Methods useful for manipulating intervals in an LAV file."""
__version__ = '0.8'

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


def conversion_coefficients(block, query_to_subject=True):
  """Return the coefficients needed to convert a query coordinate in "block"
  into a subject coordinate.
  Returns a tuple of (strand, offset), where
  subject_coord = query_coord * strand + offset
  strand is -1 if the query is reverse complemented, 1 otherwise
  """
  origin = 'query'
  target = 'subject'
  if not query_to_subject:
    (origin, target) = (target, origin)
  hit_origin = getattr(block.parent.parent, origin)
  if hit_origin['revcomp']:
    strand = -1
  else:
    strand = 1
  block_target = getattr(block, target)
  block_origin = getattr(block, origin)
  offset = block_target['begin'] - strand * block_origin['begin']
  return (strand, offset)


def hits_to_conv_table(hits, contigs=None, tuples=False, query_to_subject=True):
  """Build a coordinate conversion table from an LAV file.
  The table lists the values needed to convert sites in each gap-free block
  from query coordinates to subject coordinates.
  "hits" can be a list of LavHits or a full LavReader object.
  "contigs" lists the names of the valid query sequences to add to the table.
    If "contigs" is given, any hit whose query name is not in "contigs" will be
    left out of the table.
  The return value is a list of dicts, one per LavBlock. The key/values of each
  dict are:
    chrom1: the chromosome the block is in (the query name)
    begin:  the block's start coordinate (in the query)
    end:    the block's end coordinate (in the query)
    chrom2: the name of the corresponding subject chromosome
    strand: the strand value for coordinate conversion
    offset: the offset for coordinate conversion
    id:     the % identity of the block
    score:  the score of the block's parent alignment
  "begin" will always be <= "end".
  "strand" and "offset" are the output of conversion_coefficients().
  If "tuples" is True, it will return a list of tuples instead, with the
  elements in the same order as above.
  If "query_to_subject" is False, it build a table for converting in the
  opposite direction (swap "query" with "subject" above).
  """
  origin = 'query'
  target = 'subject'
  if not query_to_subject:
    (origin, target) = (target, origin)
  table = []
  for hit in hits:
    chrom1 = getattr(hit, origin)['name'].split()[0]
    if contigs is not None and chrom1 not in contigs:
      continue
    new_blocks = alignments_to_conv_table(hit,
                                          query_to_subject=query_to_subject)
    if tuples:
      for b in new_blocks:
        block_tuple = (b['chrom1'], b['begin'], b['end'], b['chrom2'],
                       b['strand'], b['offset'], b['id'], b['score'])
        table.append(block_tuple)
    else:
      table.extend(new_blocks)
  return table


def alignments_to_conv_table(alignments, query_to_subject=True):
  origin = 'query'
  target = 'subject'
  if not query_to_subject:
    (origin, target) = (target, origin)
  table = []
  chrom1 = chrom2 = None
  for alignment in alignments:
    if chrom1 is None:
      chrom1 = getattr(alignment.parent, origin)['name'].split()[0]
    if chrom2 is None:
      chrom2 = getattr(alignment.parent, target)['name'].split()[0]
    for block in alignment:
      block_origin = getattr(block, origin)
      begin = block_origin['begin']
      end = block_origin['end']
      if begin > end:
        (begin, end) = (end, begin)
      (strand, offset) = conversion_coefficients(block,
                                              query_to_subject=query_to_subject)
      id_ = block.identity
      score = block.parent.score
      table.append({'chrom1':chrom1, 'begin':begin, 'end':end, 'chrom2':chrom2,
                    'strand':strand, 'offset':offset, 'id':id_, 'score':score})
  return table


def convert(table, coord, chrom=None, fail='giveup', choose='id'):
  """Convert a coordinate using a conversion table produced by
  blocks_to_conv_table().
  Current implementation is naive: Compares the coordinate to every block in
  the table to check if it's contained in it.
  If "chrom" is not given, it will use any block whose interval includes
  "coord", regardless of the chromosome name."""
  assert fail in ('throw', 'giveup', 'tryharder'), 'Invalid "fail" parameter.'
  assert choose in ('id', 'score', 'length'), 'Invalid "choose" parameter.'
  assert len(table) == 0 or isinstance(table[0], dict), (
    '"table" must contain dicts, not tuples.')
  # Get a list of all blocks that contain the site.
  containing_blocks = []
  for block in table:
    if ((chrom is None or chrom == block['chrom1']) and
        block['begin'] <= coord <= block['end']):
      containing_blocks.append(block)
  # Get best_block: the best choice given the outcome of the above search.
  # matched no blocks
  if len(containing_blocks) == 0:
    if fail == 'throw':
      raise Exception('No matching blocks for the given coordinate.', 'fail')
    elif fail == 'giveup':
      return (coord, chrom)
    elif fail == 'tryharder':
      (left_block, right_block) = _flanking_blocks(table, coord, chrom)
      if left_block is None and right_block is None:
        return (coord, chrom)
      else:
        return _interpolate_coord(coord, left_block, right_block)
  # matched one block
  elif len(containing_blocks) == 1:
    best_block = containing_blocks[0]
  # matched multiple blocks
  else:
    best_block = _choose_block(containing_blocks, choose)
  # Do the actual conversion
  return convert_with_block(coord, best_block)


def convert_with_block(coord, block):
  """The actual coordinate conversion calculation, once a block is selected."""
  return (coord * block['strand'] + block['offset'], block['chrom2'])


def _choose_block(blocks, choose):
  """Choose between blocks according to the given criteria."""
  best_block = blocks[0]
  for block in blocks:
    if choose == 'id':
      if block['id'] > best_block['id']:
        best_block = block
    elif choose == 'score':
      if block['score'] > best_block['score']:
        best_block = block
    elif choose == 'length':
      if block['end'] - block['begin'] > best_block['end'] - best_block['begin']:
        best_block = block
    else:
      raise AssertionError('"choose" must be "id", "score", or "length"')
  return best_block


def _flanking_blocks(table, coord, chrom):
  """Find the closest blocks on either side of the given coordinate.
  This ignores blocks that actually contain the coordinate.
  If there are no blocks on one side of the coordinate, None will be returned
  for that side.
  If two blocks are equally close, it chooses the first one in the table.
  Can still fail if "chrom" does not match the chrom of any of the blocks.
  Will return None in that case."""
  left_best_block = None
  left_best_distance = 4294967295
  right_best_block = None
  right_best_distance = 4294967295
  for block in table:
    if chrom is None or chrom == block['chrom1']:
      # block is left of coord
      if block['end'] < coord:
        left_distance = coord - block['end']
        if left_distance < left_best_distance:
          left_best_distance = left_distance
          left_best_block = block
      # block is right of coord
      elif coord < block['begin']:
        right_distance = block['begin'] - coord
        if right_distance < right_best_distance:
          right_best_distance = right_distance
          right_best_block = block
  return (left_best_block, right_best_block)


def _interpolate_coord(coord, left_block, right_block):
  """Convert the coordinate of a site that doesn't exist in an aligned region.
  This could either be because of an insertion in the origin sequence or a
  region in the origin that just doesn't align to the target.
  This will pick a coordinate in-between the nearest aligned bases in the
  target, proportional to the distance from each.
  If one block is None, it will use the other block to convert, as if it
  extended out to the coordinate.
  This will return an answer even if the blocks are in opposite orientations
  (the coordinate of the closest base in the closest block)."""
  assert left_block is not None or right_block is not None, (
    'One of the blocks must not be None.')
  # If one of the blocks is missing, just convert according to the other block.
  if left_block is None:
    return convert_with_block(coord, right_block)
  elif right_block is None:
    return convert_with_block(coord, left_block)
  # Both blocks are present
  assert left_block['end'] < right_block['begin'], (
    'The left block must be to the left of the right block!')
  new_chrom = left_block['chrom2']
  # If orientation of the blocks don't match, use the coord of the closest base.
  if left_block['strand'] != right_block['strand']:
    # the coordinate is closer to the left block than the right block
    if coord - left_block['end'] < right_block['begin'] - coord:
      new_coord = (left_block['end'] * left_block['strand'] +
                   left_block['offset'])
    else:
      new_coord = (right_block['end'] * right_block['strand'] +
                   right_block['offset'])
  # Normal case: both blocks are present and in the same orientation
  else:
    loffset = left_block['offset']
    roffset = right_block['offset']
    lend = left_block['end']
    rend = right_block['begin']
    # calculate an intermediate offset between those of the two blocks
    offset = loffset + int(round((roffset - loffset) *
                                 ((coord - lend)/(rend - lend))))
    new_coord = coord * left_block['strand'] + offset
  return (new_coord, new_chrom)


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