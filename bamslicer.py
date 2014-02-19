#!/usr/bin/env python
"""Methods to select reads from a BAM that contain given variants, and return
statistics on their characteristics."""
from __future__ import division
import os
import sys
import collections
from pyBamParser.bam import Reader

NUM_FLAGS = 12
DEFAULT_MAX_MAPQ = 40
STAT_NAMES = ['supporting', 'coverage', 'flags', 'mapqs', 'freq', 'strand_bias',
  'mate_bias']

def get_reads_and_stats(bamfilepath, variants, supporting=True, opposing=False,
    readgroups=None):
  """Take a BAM and a list of variants, and return the reads covering the
  variants, as well as statistics on those reads.
  ***Note: requesting both supporting and opposing reads is not allowed at the
  moment.***
  The selected reads support and/or oppose the given variants, depending on the
  value of those arguments (if both supporting and opposing are True, it will
  select all reads covering that site). The statistics summarize whatever reads
  are selected.
  The BAM file and variants list should be sorted by start coordinate. Each
  variant is a dict, e.g.:
    {'chrom':'chrM', 'coord':310, 'type':'S', 'alt':None}
    {'chrom':'1', 'coord':2345, 'type':'D', 'alt':'2'}
    {'chrom':'pUC18', 'coord':4210, 'type':'I', 'alt':'GAT'}
  Return value: A 2-tuple of reads and stats. Both are lists where each element
  corresponds to a variant of the same index in the input variants list. Each
  element in the reads list is a list of BAMReads which support and/or oppose
  that variant. The elements of the stats list summarize each set of reads. See
  the description of get_read_stats() for the format.
  """
  #TODO: see about chromosome issues in order of reads and variants
  if supporting and opposing:
    raise NotImplementedError

  bam_reader = Reader(bamfilepath)

  read_sets_supporting = [None] * len(variants)
  read_sets_opposing = [None] * len(variants)
  stat_sets = [None] * len(variants)
  for i in range(len(variants)):
    read_sets_supporting[i] = []
    read_sets_opposing[i] = []
    stat_sets[i] = {}

  for read in bam_reader:
    read_rname  = read.get_reference_name()
    read_pos    = read.get_position()
    read_end    = read.get_end_position()
    (read_ins, read_del) = read.get_indels()

    #TODO: use deque for local variants by discarding variants we've moved past
    #      and only adding ones to the end when needed
    #current_variants = collections.deque()
    for (i, variant) in enumerate(variants):
      var_pos = variant['coord']
      if not (read_pos <= var_pos <= read_end):
        continue
      var_chrom = variant['chrom']
      if read_rname != var_chrom:
        continue
      # now check if the read supports the variant
      #TODO: add support for checking the alt allele identity
      var_type = variant['type']
      supports = False
      if var_type == 'I':
        for (ins_pos, ins_len) in read_ins:
          if ins_pos == var_pos:
            supports = True
      elif var_type == 'D':
        for (del_pos, del_len) in read_del:
          if del_pos == var_pos:
            supports = True
      else:
        raise NotImplementedError('Unsupported variant type "'+var_type+'"')
      # add read to the appropriate list of reads
      if supports:
        read_sets_supporting[i].append(read)
      else:
        read_sets_opposing[i].append(read)

  for (i, reads_supporting) in enumerate(read_sets_supporting):
    stat_sets[i] = get_read_stats(reads_supporting, variants[i],
      reads_opposite=read_sets_opposing[i])

  if supporting:
    return (read_sets_supporting, stat_sets)
  else:
    return (read_sets_opposing, stat_sets)


#TODO: Break this into individual functions that add one stat each. No more
#      reads/reads_opposite distinction, no more stats_to_get checks.
def get_read_stats(reads, variant, reads_opposite=None, stats_to_get=STAT_NAMES):
  """Calculate a set of summary statistics for the set of reads.
  If stats_to_get is given, it will only include the statistics whose keys are
  provided in that list. Other statistics will be None. reads_opposite is a list
  of reads that oppose the variant, if reads is the list of reads supporting it.
  If reads are the reads opposing the variant, then vice-versa.
  Returns a dict. Descriptions of the keys and their values:
    'supporting':
  The number of supporting reads.
    'coverage':
  Total number of reads spanning this position.
  None if no opposite reads are provided.
    'flags':
  The totals of how many reads have each flag set.
  Format: a list, where flags[i] = the total number of reads that have the
  flag 2**i set.
    'mapqs':
  The totals of how many reads have each MAPQ value.
  Format: a list, where mapq[value] = the total number of reads that have a
  MAPQ equal to "value". Note the list is 0-based, so MAPQ value 40 will be the
  41st element in it.
    'freq':
  The frequency of the variant.
  freq = len(reads)/(len(reads)+len(reads_opposing))
  None if no opposite reads are provided.
    'strand_bias':
  Strand bias of the variant.
  Based on method 1 (SB) of Guo et al., 2012.
  None if no opposite reads are provided.
    'mate_bias':
  Mate bias of the variant toward 1st or 2nd read in the pair.
  Calculated in same way as strand_bias.
    'read_pos':
  The distribution of where the variant occurs along the length of the reads.
  Note: if the reads given are the opposing reads, it will still look for the
  variant in them and report its findings (zeroes), wasting time. Consider not
  asking for this statistic in that case.
    'flank_quals':
  ***PLANNED*** The PHRED quality of the bases flanking the
  variant
    'nm_edits':
  ***PLANNED*** NM tag edit distances.
  """
  #TODO: doublecheck values that could be incalculable given the input data
  stats = {}

  if 'supporting' in stats_to_get:
    stats['supporting'] = len(reads)
  else:
    stats['supporting'] = None

  if 'coverage' in stats_to_get and reads_opposite is not None:
    stats['coverage'] = len(reads) + len(reads_opposite)
  else:
    stats['coverage'] = None

  if 'mapqs' in stats_to_get:
    stats['mapqs'] = sum_mapqs(reads)
  else:
    stats['mapqs'] = None

  if 'flags' in stats_to_get:
    stats['flags'] = sum_flags(reads)
  else:
    stats['flags'] = None

  if 'strand_bias' in stats_to_get and reads_opposite is not None:
    if stats['flags'] is None:
      flags_for = sum_flags(reads)
    else:
      flags_for = stats['flags']
    flags_against = sum_flags(reads_opposite)
    stats['strand_bias'] = get_strand_bias(flags_for, flags_against,
      len(reads), len(reads_opposite))
  else:
    stats['strand_bias'] = None

  if 'mate_bias' in stats_to_get and reads_opposite is not None:
    if stats['flags'] is None:
      flags_for = sum_flags(reads)
    else:
      flags_for = stats['flags']
    flags_against = sum_flags(reads_opposite)
    stats['mate_bias'] = get_mate_bias(flags_for, flags_against)
  else:
    stats['strand_bias'] = None

  if ('freq' in stats_to_get and
      reads_opposite is not None and
      len(reads) + len(reads_opposite) > 0):
    stats['freq'] = len(reads)/(len(reads)+len(reads_opposite))
  else:
    stats['freq'] = None

  if 'read_pos' in stats_to_get:
    stats['read_pos'] = get_read_pos(reads, variant)

  return stats


def sum_mapqs(reads):
  mapqs = [0] * (DEFAULT_MAX_MAPQ + 1)
  for read in reads:
    read_mapq = read.get_mapq()
    if read_mapq > len(mapqs) - 1:
      mapqs.extend([0] * (read_mapq - len(mapqs) + 1))
    mapqs[read_mapq]+=1
  return mapqs


def sum_flags(reads):
  flags_sum = [0] * NUM_FLAGS
  for read in reads:
    set_flags = get_flags(read.get_flag())
    for i in range(NUM_FLAGS):
      if set_flags[i]:
        flags_sum[i]+=1
  return flags_sum


def get_flags(flagint):
  """Give a SAM flag as an integer, return a list of whether each flag is set.
  Returns a list of booleans, where set_flags[i] = whether flag 2**i is set"""
  flagbin = ('{0:0'+str(NUM_FLAGS)+'b}').format(flagint)
  i = 0
  set_flags = [False] * NUM_FLAGS
  for bit in reversed(flagbin):
    set_flags[i] = bit == '1'
    i+=1
  return set_flags


def get_strand_bias(flags_for, flags_against, total_for, total_against):
  """Based on method 1 (SB) of Guo et al., 2012.
  If there a denominator that would be 0, there is no valid result and this will
  return None. This occurs when there are no reads on one of the strands, or
  when there are no minor allele reads."""
  a = total_against - flags_against[4]
  b = total_for - flags_for[4]
  c = flags_against[4]
  d = flags_for[4]
  return get_bias(a, b, c, d)


def get_mate_bias(flags_for, flags_against):
  """Based on method 1 (SB) of Guo et al., 2012.
  If there a denominator that would be 0, there is no valid result and this will
  return None. This occurs when there are no reads from one of the mates, or
  when there are no minor allele reads."""
  a = flags_against[6]
  b = flags_for[6]
  c = flags_against[7]
  d = flags_for[7]
  return get_bias(a, b, c, d)


def get_bias(a, b, c, d):
  """Bare equation from Guo et al., 2012 (SB, strand bias method 1).
  One modification: a/b and c/d will be swapped if necessary to make sure b and
  d are the minor allele..
  """
  if b + d > a + c:
    (a, b, c, d) = (b, a, d, c)
  try:
    return abs(b/(a+b) - d/(c+d)) / ((b+d) / (a+b+c+d))
  except ZeroDivisionError:
    return None


def get_read_pos(reads, variant):
  """Get the variant's distribution of positions along the reads.
  Returns a dict of the offset of the variant in each read. Keys are the offsets
  (variant coord - read start coord) and values are counts (how many reads had
  the variant at that position)."""
  raise NotImplementedError
  read_pos = {}
  for read in reads:
    pass
  return read_pos