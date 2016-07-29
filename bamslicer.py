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

def get_variant_stats(bamfilepath, variants, ref=None):
  """Take a BAM and a dict of lists of variants, and return statistics on the
  reads covering each.
  The BAM file should be sorted by chromosome, then coordinate (the standard
  "samtools sort"), and the variants should be a dict mapping chromosome names to
  lists of the variants on each chromosome, sorted by start coordinate.
  Each variant is a dict, e.g.:
    {'chrom':'1', 'coord':2345, 'type':'D', 'alt':'2'}
    {'chrom':'pUC18', 'coord':4210, 'type':'I', 'alt':'GAT'}
    {'chrom':'chrM', 'coord':310, 'type':'S', 'alt':None}
  NOTE: Currently SNV's are not supported!
  Return value: A list of read statistics.
  Each value in the list represents the statistics for one variant.
  Each value is a dict mapping sample (read group) names to statistics for the
  variant in that sample. Each set of statistics is a dict mapping the names of
  different statistics to their values. For convenience, the dicts include the
  original variant information (the 'chrom', 'coord', etc keys above).
  In addition, they currently include these statistics keys:
  "supporting", "opposing", "mapqs", "flags", "flags_opposing", "strand_bias",
  "mate_bias", "var_pos_dist", and "context" when a reference is given.
  N.B.: When inspecting reads, this cannot distinguish between insertions of the
  same length, but different inserted sequence. So this could falsely identify
  reads as supporting a particular insertion when in fact they support a different
  insertion of the same length. This is due to a limitation of pyBamParser (it
  currently doesn't correlate the CIGAR string field with the sequence field).
  """
  # This function relies on both the input BAM and list of variants being sorted by coordinate.
  # It reads through this sorted BAM one read at a time. As it reads, it keeps a list of variants
  # inside the footprint of the current read. They way it does this is that as it takes on each new
  # read, it drops variants that are now to the left of the current read, and adds new variants
  # which are now covered by the right side of the read. Then it checks each of these variants
  # against each of the indels in the read.

  if 'TEST_OUTPUT' in os.environ:
    test_file = open(os.path.expanduser(os.environ['TEST_OUTPUT']), 'w')
  else:
    test_file = None

  bam_reader = pyBamParser.bam.Reader(bamfilepath)

  last_chrom = None
  for read in bam_reader:

    #TODO: Check that bam is sorted by chromosome, then coordinate.
    read_chrom = read.get_reference_name()
    read_start = read.get_position()
    read_end   = read.get_end_position()

    # If we're starting a new chromosome, re-initialize all variables.
    if read_chrom != last_chrom:
      chrom_variants = variants.get(read_chrom, [])
      current_variants = []
      reads_supporting = collections.defaultdict(dict)
      reads_opposing = collections.defaultdict(dict)
      current_start = read_start
      current_end = read_start
      variants_i = 0
      last_chrom = read_chrom

    # Pop variants off the left end and process them.
    while (current_start < read_start and
           current_variants and current_variants[0]['coord'] < read_start):
      variant = current_variants.pop(0)
      key = (variant['coord'], variant['type'], variant['alt'])
      if test_file:
        _test_write(test_file, variant, reads_opposing[key])
      yield get_samples_read_stats(variant, reads_supporting[key], reads_opposing[key], ref)
      del(reads_supporting[key])
      del(reads_opposing[key])
      if current_variants:
        current_start = current_variants[0]['coord']
        #TODO: Something smarter when the list is empty? This could happen when there's a big gap
        #      between reads and we reach the end of the last one before the start of this one.

    # Add new variants to the right end.
    if current_end < read_end:
      while variants_i < len(chrom_variants) and chrom_variants[variants_i]['coord'] <= read_end:
        variant = chrom_variants[variants_i]
        #TODO: We might skip some variants here (if they fall in gaps between reads). Should report
        #      something about these variants? Basically that there was no information?
        if variant['coord'] >= read_start:
          current_variants.append(variant)
          current_end = variant['coord']
        variants_i += 1

    # Make sure current_start is the first variant.
    # Or, if there are none in this read, skip the rest.
    if current_variants:
      current_start = current_variants[0]['coord']
    else:
      continue

    # Get sample name from read group (returns None if no RG tag).
    sample = read.get_read_group()
    (read_ins, read_del) = read.get_indels()

    for variant in current_variants:
      var_pos = variant['coord']
      # Double-check that the variant is in the footprint of the read.
      if not (read_start <= var_pos <= read_end):
        continue
      # Does the read support the variant?
      var_type = variant['type']
      var_alt = variant['alt']
      supports = False
      if var_type == 'I':
        for ins_pos, ins_len in read_ins:
          if ins_pos == var_pos:
            if var_alt is None or ins_len == len(var_alt):
              supports = True
      elif var_type == 'D':
        for del_pos, del_len in read_del:
          if del_pos == var_pos:
            if var_alt is None or del_len == int(var_alt):
              supports = True
      else:
        raise NotImplementedError('Unsupported variant type "'+var_type+'"')
      # Add read to the appropriate list of reads.
      # If sample not in dict for this variant, add it.
      key = (variant['coord'], variant['type'], variant['alt'])
      if supports:
        #TODO: Figure out how to make second-tier dict a defaultdict?
        reads_supporting_list = reads_supporting[key].get(sample, [])
        reads_supporting_list.append(read)
        reads_supporting[key][sample] = reads_supporting_list
      else:
        #TODO: Figure out how to make second-tier dict a defaultdict?
        reads_opposing_list = reads_opposing[key].get(sample, [])
        reads_opposing_list.append(read)
        reads_opposing[key][sample] = reads_opposing_list

  # Process the remaining variants.
  for variant in current_variants:
    key = (variant['coord'], variant['type'], variant['alt'])
    if test_file:
      _test_write(test_file, variant, reads_opposing[key])
    yield get_samples_read_stats(variant, reads_supporting[key], reads_opposing[key], ref)


def get_samples_read_stats(variant, sample_reads_supporting, sample_reads_opposing, ref=None):
  """Compute read statistics for the variant in each sample.
  Return a dict mapping sample names to statistics dicts."""
  sample_stats = {}
  for sample in sample_reads_supporting:
    reads_supporting = sample_reads_supporting[sample]
    reads_opposing   = sample_reads_opposing[sample]
    var_stats = get_read_stats(reads_supporting, reads_opposing)
    var_stats['var_pos_dist'] = get_var_pos_dist(reads_supporting, variant)
    if ref:
      var_stats['context'] = get_context(variant, ref)
    var_stats.update(variant)
    sample_stats[sample] = var_stats
  return sample_stats


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
  stats = {}

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


#TODO: use samflags.py
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


#TODO: Optimize. Currently, for every variant, this opens the ref file, reads through the entire
#      reference until the variant's position, and closes the ref file again.
def get_context(variant, refpath, flank_len=DEFAULT_FLANK_LEN):
  """Return the (reference) sequence context surrounding the variant.
  The returned sequence will be 2*flank_len long*, centered on variant['coord'].
  E.g. for {'chrom':'chrM','coord':12417,'type':'S','alt':'G'} it will give
  'CTCGTTAACCCTAAcAAAAAAAACTCATAC', where the lowercased "c" is base 12417,
  with 14bp of left flank sequence and 15bp of right flank.
  If the variant is an insertion, the base after variant['coord'] will be lower-
  cased. The same goes for a deletion with a None variant['alt']. If an alt is
  provided, the deleted sequence will be lowercased.
  *If the edge of the chromosome is closer than flank_len, the output will be
  shorter, consisting of the bases up to the end."""
  if refpath is None:
    return None
  assert variant['type'] in 'SID', 'Variant type must be one of "S", "I", "D".'
  # get the sequence flank_len bp before and after the variant
  start = max(1, variant['coord'] - flank_len + 1)
  end = variant['coord'] + flank_len
  fasta = fastareader.FastaLineGenerator(refpath)
  # print "requesting %s:%s-%s" % (variant['chrom'], start, end)
  seq = fasta.extract(start, end, chrom=variant['chrom'])
  # lowercase the variant
  lflank_len = variant['coord'] - start + 1
  if variant['type'] in 'S':
    var_start = lflank_len - 1
    var_end = var_start + 1
  elif variant['type'] == 'D':
    var_start = lflank_len
    if variant['alt'] is None:
      var_end = var_start + 1
    else:
      var_end = var_start + int(variant['alt'])
  elif variant['type'] == 'I':
    var_start = lflank_len
    var_end = var_start + 1
  lflank = seq[:var_start]
  var_seq = seq[var_start:var_end].lower()
  rflank = seq[var_end:]
  return lflank+var_seq+rflank


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


def _test_write(test_file, variant, reads_opposing):
  all_reads = []
  for reads_list in reads_opposing.values():
    for read in reads_list:
      all_reads.append(read)
  reads = ','.join([read.get_read_name() for read in all_reads])
  test_file.write('{coord}\t{type}\t{alt}\t{reads}\n'
                   .format(reads=reads, **variant))

version_check(EXPECTED_VERSIONS)
