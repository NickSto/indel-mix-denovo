#!/usr/bin/env python
import sys
import pysam

DEFAULT_CHROM = 'chrM'
DEFAULT_SIZE = 16569
PHRED_THRES = 30
MAPQ_THRES = 20
ALLOWED_FLAGS = (147,99,163,83,1107,1187,1123,1171)

scriptname = sys.argv[0].split('/')[-1]
USAGE = 'Usage: $ '+scriptname+""" alignment.bam
Consensus sequence will be written to alignment.bam.major.fa
"""

chrom = DEFAULT_CHROM
size = DEFAULT_SIZE

# Get SAM/BAM file from command line
if len(sys.argv) > 1:
  sam = sys.argv[1]
else:
  sys.stderr.write(USAGE)
  sys.exit(1)

fasta = sam+'.major.fa'
samfile = pysam.Samfile(sam, 'rb')

def select_read(read):
  """Return the read if it passes the thresholds, None otherwise."""
  if (read.alignment.mapq >= MAPQ_THRES and
      read.alignment.flag in ALLOWED_FLAGS):
    return read
  else:
    pass

def get_major(bases):
  """Return the most common base in the list of bases.
  Input is a list of single characters, 'A', 'C', 'G', or 'T'."""
  major = {0:'A',1:'C',2:'G',3:'T'}
  counts = [
    bases.count('A'),
    bases.count('C'),
    bases.count('G'),
    bases.count('T')
  ]
  return major[counts.index(max(counts))] 

def printseq(seq, outfile, seqname):
  """Print out sequence to FASTA file"""
  print >> outfile, '>'+seqname
  for i in range(0, len(seq), 70):
    print >> outfile, ''.join(seq[i:i+70])

sequence = ['N'] * size

# For each position in the reference:
for pileupcolumn in samfile.pileup(chrom, 0, size, stepper='all', max_depth=10000000, mask=False):
  pos_in_ref = pileupcolumn.pos
  # The list bases in each read at this position (a list of single characters)
  bases_in_position = []
  # For each read at this position:
  for pileupread in pileupcolumn.pileups:
    # Does the read pass MAPQ and flags thresholds?
    if select_read(pileupread):
      alignment = pileupread.alignment
      pos_in_read = pileupread.qpos
      # Does the base pass the PHRED quality threshold?
      if ord(alignment.qual[pos_in_read])-33 >= PHRED_THRES:
        # Then count it in the list of bases at this position
        bases_in_position.append(alignment.seq[pos_in_read])
  # Find the major allele and add it to the sequence
  sequence[pos_in_ref] = get_major(bases_in_position)

# Write consensus sequence to output FASTA file
with open(fasta, 'w+') as fastafile:
  printseq(sequence, fastafile, chrom)
