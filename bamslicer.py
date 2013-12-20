#!/usr/bin/env python
"""Methods to select reads from a BAM that contain given variants, and return
statistics on their characteristics."""
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
import collections
import sys
import os

NUM_FLAGS = 12
DEFAULT_MAX_MAPQ = 40
STAT_NAMES = ['flags', 'mapqs']

"""End uses of this library:

nvc-filter.py:
Give a variant and a bam file (or a bam_reader), and return statistics on the
reads containing (and not containing?) the variant.
If it could do this without nvc-filter.py having to import any pyBamParser
stuff, that'd be neat.

inspect-reads.py:
Give a variant and a bam_reader, and return the reads containing, not
containing, or covering the variants, and statistics on those reads.
"""

"""Design:

This is a problem that really needs to take advantage of BAM indexing. The main
issue is figuring out which reads overlap with each variant, and that's what
the index's bins are made for.

So for now, a totally naive implementation. It will have to go through the
entire BAM file on each call. This makes it basically impossible to use by
nvc-filter.py, which has to call it once for every variant it investigates.
"""


def get_reads(bamfilepath, variants, supporting=True, opposing=False):
  """Take a BAM and a list of variants, and return the reads covering the
  variants. The BAM and variants list should be sorted by start coordinate."""
  # The sorting of the BAM and variants list is crucial to keeping this from
  # taking fully quadratic time. It still is, really, but this way I effectively
  # only look at variants between the start and end of the current read.
  #TODO: see about chromosome issues in order of reads and variants

  bam_reader = Reader(bamfilepath)
  

  reads = []
  for read in bam_reader:
    read_pos    = read.get_position()
    read_end    = read.get_end_position()
    read_rname  = read.get_reference_name()
    read_indels = read.get_indels()

    #TODO: use deque for local variants by discarding variants we've moved past
    #      and only adding ones to the end when needed
    #current_variants = collections.deque()
    for variant in variants:
      var_pos = variant['coord']
      if not (read_pos <= var_pos <= read_end):
        continue
      var_chrom = variant['chrom']
      if read_rname != var_chrom:
        continue
      if supporting and opposing:
        reads.append(read)
        break
      # now check if the read supports the variant
      var_type = variant['type']
      supports = False
      if var_type == 'I':
        for (ins_pos, ins_len) in read_indels[0]:
          if ins_pos == var_pos:
            supports = True
      elif var_type == 'D':
        for (del_pos, del_len) in read_indels[1]:
          if del_pos == var_pos:
            supports = True
      else:
        raise Exception('Unsupported variant type "'+var_type+'"')
      if (supporting and supports) or (opposing and not supports):
        reads.append(read)
        break

  return reads


def get_read_stats(reads, variant, stats_to_get=STAT_NAMES):
  """Calculate a set of summary statistics for the set of reads.
  Returns a dict. Descriptions of the keys and their values:
    'flags':
  The totals of how many reads have each flag set.
  Format: a list, where flags[i] = the total number of reads that have the
  flag 2**i set.
    'mapqs':
  The totals of how many reads have each MAPQ value.
  Format: a list, where mapq[value] = the total number of reads that have a
  MAPQ equal to "value". Note the list is 0-based, so MAPQ value 40 will be the
  41st element in it.
    'read_pos':
  ***PLANNED*** The distribution of where the variant occurs along the length of
  the reads
    'flank_quals':
  ***PLANNED*** The PHRED quality of the bases flanking the
  variant
    'nm_edits':
  ***PLANNED*** NM tag edit distances.
  """

  stats = {}
  stats['flags'] = [0] * NUM_FLAGS
  stats['mapqs'] = [0] * (DEFAULT_MAX_MAPQ + 1)
  for read in reads:

    if 'flags' in stats_to_get:
      set_flags = get_set_flags(read.get_flag())
      for i in range(NUM_FLAGS):
        if set_flags[i]:
          stats['flags'][i]+=1

    if 'mapqs' in stats_to_get:
      read_mapq = read.get_mapq()
      if read_mapq > len(stats['mapqs']) - 1:
        stats['mapqs'].extend([0] * (read_mapq - len(stats['mapqs']) + 1))
      stats['mapqs'][read_mapq]+=1

  return stats


def get_set_flags(flagint):
  flagbin = ('{0:0'+str(NUM_FLAGS)+'b}').format(flagint)
  i = 0
  set_flags = [False] * NUM_FLAGS
  for bit in reversed(flagbin):
    set_flags[i] = bit == '1'
    i+=1
  return set_flags