#!/usr/bin/env python
"""Methods to select reads from a BAM that contain given variants, and return
statistics on their characteristics."""
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
import collections
import sys
import os


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

Both uses will be inspecting one variant at a time, but in sequential order, in
the same BAM file. Going through all the reads every time is a non-starter,
especially for nvc-filter.py.

How about making an object that you add variants to, and it moves through the
BAM file as you add them, discarding reads before the current variants?
"""

class Variant(object):
  """Simple representation of a variant, mainly to make it hashable.
  Types for attributes: chrom = str, pos = int, type = str, alt = str
  type should be one of 'S', 'I', or 'D' for SNV, insertion, or deletion.
  alt is optional, and should be the alternate base for SNVs, the inserted
  sequence for insertions (INCLUDING the preceding base, as in VCF), or the
  length of the deletion."""

  def __init__(self, chrom, pos, vartype, alt=None):
    self.chrom = chrom
    self.pos = pos
    self.type = vartype
    self.alt = alt

  # representation as a tuple, for equality and hashing purposes
  def __key(self):
    return (self.chrom, self.pos, self.type, self.alt)

  def __eq__(a, b):
    return type(a) == type(b) and a.__key() == b.__key()

  def __hash__(self):
    return hash(self.__key())


class BamVariants(object):
  """Go through a list of Variants, building information on their presence in
  a BAM file. The BAM file must be sorted, and the Variants must be added in
  order (of their pos)."""

  def __init__(self, bamfilepath):
    #TODO: a check on whether the BAM file is sorted
    self._bam_reader = Reader(bamfilepath)
    self._read_buffer = collections.deque()
    self._read_buffer_ends = collections.deque()
    self._next_read
    self._stats = {}
    self._reads = {}
    self._last_stats = {}
    self._last_reads = []


  def add_variant(self, variant):
    """Add a variant
    Variant format examples:
    {'chrom':'chrM', 'coord':310, 'type':'S', 'alt':None}
    {'chrom':'1', 'coord':2345, 'type':'D', 'alt':'2'}
    {'chrom':'pUC18', 'coord':4210, 'type':'I', 'alt':'GAT'}
    """

    var_pos = variant.pos
    done = False
    while self._read_buffer:
      read_end = self._read_buffer_ends[0]
      if variant.pos < read.get_position():
        # we're past this 
        self._read_buffer.popleft()


  def get_stats(self, variant):
    return self._stats.get(variant, None)

  def get_reads(self, variant):
    return self._reads.get(variant, None)

  def get_last_stats(self):
    """Return the statistics of the last variant added"""
    return self._last_stats

  def get_last_reads(self, supporting=True, opposing=False):
    """Return the reads associated with the last variant added"""
    return self._last_reads




def get_reads(bamfilepath, variants, supporting=True, opposing=False):
  """Take a BAM and a list of variants, and return the reads covering the
  variants. The BAM and variants list should be sorted by start coordinate."""
  # The sorting of the BAM and variants list is crucial to keeping this from
  # taking fully quadratic time. It still is, really, but this way I effectively
  # only look at variants between the start and end of the current read.
  #TODO: see about chromosome issues in order reads and variants

  bam_reader = Reader(bamfilepath)
  variants_deque = collections.deque(variants)

  reads = []
  for read in bam_reader:
    read_pos   = read.get_position()
    read_end   = read.get_end_position()
    read_rname = read.get_reference_name()
    # read_qname = read.get_read_name()
    # read_cigar = read.get_sam_cigar()

    popleft = 0
    for variant in variants_deque:
      var_pos = variant['coord']
      if var_pos < read_pos:
        # drop this variant from the queue: we're at reads beyond it
        popleft += 1
      elif var_pos > read_end:
        # finished: we're at variants beyond the current read
        break
      var_chrom = variant['chrom']
      if read_rname != var_chrom:
        continue
      # now check if the read contains the variants that are left







