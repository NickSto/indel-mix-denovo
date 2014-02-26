#!/usr/bin/env python
"""Methods to select reads from a BAM that contain given variants, and return
statistics on their characteristics."""
from __future__ import division
import os
import sys
import collections
import pkg_resources
import pyBamParser.bam
import fastareader
__version__ = '0.6'

EXPECTED_VERSIONS = {'fastareader':'0.7', 'pyBamParser':'0.0.1'}

NUM_FLAGS = 12
DEFAULT_MAX_MAPQ = 40
DEFAULT_FLANK_LEN = 15

def get_reads_and_stats(bamfilepath, variants, ref=None, readgroups=None):
  """Take a BAM and a list of variants, and return the reads supporting the
  variants, as well as statistics on the reads covering it.
  The BAM file and variants list should be sorted by start coordinate. Each
  variant is a dict, e.g.:
    {'chrom':'1', 'coord':2345, 'type':'D', 'alt':'2'}
    {'chrom':'pUC18', 'coord':4210, 'type':'I', 'alt':'GAT'}
    {'chrom':'chrM', 'coord':310, 'type':'S', 'alt':None}
  NOTE: Currently SNV's are not supported!
  Return value: A 2-tuple of reads and stats. Both are lists where each element
  corresponds to a variant of the same index in the input variants list. Each
  element in the reads list is a list of BAMReads which support that variant.
  The elements of the stats list summarize each set of reads. See the
  description of get_read_stats() for the format.
  N.B.: When looking for the variant in a read, this checks for the same ALT
  allele, but only the length, not the identity. This only matters for
  insertions, not deletions. So it will fail to distinguish insertions of the
  same length but different sequences.
  """

  bam_reader = pyBamParser.bam.Reader(bamfilepath)

  # pre-allocate lists of reads so they can be randomly accessed later
  read_sets_supporting = [None] * len(variants)
  read_sets_opposing = [None] * len(variants)
  for i in range(len(variants)):
    read_sets_supporting[i] = []
    read_sets_opposing[i] = []

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
      var_alt = variant['alt']
      supports = False
      if var_type == 'I':
        for (ins_pos, ins_len) in read_ins:
          if ins_pos == var_pos:
            if var_alt is None or ins_len == len(var_alt):
              supports = True
      elif var_type == 'D':
        for (del_pos, del_len) in read_del:
          if del_pos == var_pos:
            if var_alt is None or del_len == int(var_alt):
              supports = True
      else:
        raise NotImplementedError('Unsupported variant type "'+var_type+'"')
      # add read to the appropriate list of reads
      if supports:
        read_sets_supporting[i].append(read)
      else:
        read_sets_opposing[i].append(read)

  # calculate sets of stats on the reads covering each variant
  stat_sets = []
  for (variant, reads_supporting, reads_opposing) in zip(
      variants, read_sets_supporting, read_sets_opposing):
    stats = get_read_stats(reads_supporting, reads_opposing)
    stats['var_pos_dist'] = get_var_pos_dist(reads_supporting, variant)
    stats['seq'] = get_seq(variant, ref)
    stat_sets.append(stats)

  return (read_sets_supporting, stat_sets)


def get_read_stats(reads_supporting, reads_opposing):
  """Calculate a set of summary statistics for the set of reads.
  This bundles all the stats that can be obtained just from inspecting the
  supporting and opposing reads. Supplementary stats using other inputs can be
  obtained from other methods.
  'reads_supporting' is a list of reads supporting the variant.
  'reads_opposing' is a list of reads opposing the variant.
  'variant' is the variant in question, in the format described in
  get_reads_and_stats(). If not given, the 'read_pos' stat will be omitted.
  If reads are the reads opposing the variant, then vice-versa.
  Returns a dict. Descriptions of the keys and their values:
    'supporting':
  The number of supporting reads.
    'opposing':
  The number of opposing reads.
    'flags':
  The totals of how many supporting reads have each flag set.
  Format: a list, where flags[i] = the total number of reads that have the
  flag 2**i set.
    'mapqs':
  The totals of how many reads have each MAPQ value.
  Format: a list, where mapq[value] = the total number of reads that have a
  MAPQ equal to "value". Note the list is 0-based, so MAPQ value 40 will be the
  41st element in it.
    'strand_bias':
  Strand bias of the variant.
  Based on method 1 (SB) of Guo et al., 2012.
  None if no opposite reads are provided.
    'mate_bias':
  Mate bias of the variant toward 1st or 2nd read in the pair.
  Calculated in same way as strand_bias.
    'flank_quals':
  ***PLANNED*** The PHRED quality of the bases flanking the variant
    'nm_edits':
  ***PLANNED*** NM tag edit distances.
  """
  #TODO: doublecheck values that could be incalculable given the input data
  stats = collections.defaultdict(lambda: None)

  stats['supporting'] = len(reads_supporting)
  stats['opposing'] = len(reads_opposing)
  stats['mapqs'] = sum_mapqs(reads_supporting)
  stats['flags'] = sum_flags(reads_supporting)
  stats['flags_opposing'] = sum_flags(reads_opposing)
  stats['strand_bias'] = get_strand_bias(
    stats['flags'], stats['flags_opposing'],
    stats['supporting'], stats['opposing']
  )
  stats['mate_bias'] = get_mate_bias(stats['flags'], stats['flags_opposing'])
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


def get_var_pos_dist(reads, variant):
  """Get the variant's distribution of positions along the reads.
  Returns a dict of the offset of the variant in each read. Keys are the
  positions along the read, as an integer percentage of the read length (0 is
  start, 100 is end, no matter the orientation). Values are counts (how many
  reads had the variant at that position)."""
  var_pos_dist = collections.defaultdict(int)
  for read in reads:
    start = read.get_position()
    end = read.get_end_position()
    length = abs(start - end) + 1
    reverse = read.get_flag() & 0b10000 == 0b10000
    if reverse:
      var_pos = 100 * (end - variant['coord']) // length
    else:
      var_pos = 100 * (variant['coord'] - read.get_position() + 2) // length
    var_pos_dist[var_pos]+=1
  return var_pos_dist


def get_seq(variant, ref, flank_len=DEFAULT_FLANK_LEN):
  """Return the (reference) sequence context surrounding the variant."""
  if ref is None:
    return None
  chrom = variant['chrom']
  coord = variant['coord']
  


def version_check(expected):
  actual = {}
  for module_name in expected:
    try:
      actual[module_name] = pkg_resources.get_distribution(module_name).version
    except pkg_resources.DistributionNotFound:
      module = sys.modules[module_name]
      for version_name in ['version', 'VERSION', '__version__']:
        if version_name in dir(module):
          actual[module_name] = getattr(module, version_name)
  for module_name in actual:
    assert actual[module_name] == expected[module_name], (
      "Wrong version of "+module_name+". Expected: "+expected[module_name]
      +", actual: "+actual[module_name]
    )

version_check(EXPECTED_VERSIONS)
