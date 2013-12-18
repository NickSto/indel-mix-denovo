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
Give a list of variants (or a variant) and a bam_reader, and return the reads
containing, not containing, or covering the variants.
"""

"""Design:

Should probably make a bunch of low-level functions that everyone will use, and
then maybe bundle them into specialized functions that only one script or the
other will use.
"""


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
    






