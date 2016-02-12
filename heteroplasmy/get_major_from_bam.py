#!/usr/bin/env python
import sys
import pysam
import argparse
import distutils.version

# Broken by PySAM 0.8.1.
# Fails with error:
# AttributeError: 'pysam.calignmentfile.PileupRead' object has no attribute 'qpos'
# on line:
#             read.alignment.seq[read.qpos] == minor and
# Looks they probably renamed qpos to something like query_alignment_position, though they don't
# specifically mention it in the release notes:
# http://pysam.readthedocs.org/en/latest/release.html#release-0-8-1
version = pysam.__version__
if distutils.version.StrictVersion(version) >= distutils.version.StrictVersion('0.8.1'):
  raise Exception('Not compatible with PySAM versions > 0.8.0. Current version: '+version)

OPT_DEFAULTS = {'chrom':'chrM', 'length':16569, 'min_phred':30, 'min_mapq':20,
                'flags':'147,99,163,83,1107,1187,1123,1171'}

DESCRIPTION = """Create a new reference sequence based on the consensus of all the reads covering
each site. Reads and bases will be filtered out according to several thresholds: mapping quality,
PHRED score, and SAM flags. Locations where there are no bases passing the thresholds will be "N"."""

parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.set_defaults(**OPT_DEFAULTS)
parser.add_argument('bam', metavar='alignment.bam',
  help='The input BAM file.')
parser.add_argument('-o', '--output',
  help='The output FASTA file. If not given, the default is "[name of BAM].major.fa".')
parser.add_argument('-c', '--chrom', required=True,
  help='The name of the sequence to consider. Default: "%(default)s".')
parser.add_argument('-l', '--length', type=int, required=True,
  help='The length of the sequence, in bp. Default: %(default)s.')
parser.add_argument('-p', '--min-phred', type=int,
  help='The minimum PHRED score required for a base to count toward the consensus. '
       'Default: %(default)s.')
parser.add_argument('-m', '--min-mapq', type=int,
  help='The minimum MAPQ score required for a read to count. Default: %(default)s.')
parser.add_argument('-f', '--allowed-flags', dest='flags',
  help='The set of valid flags a read can have for it to count. Give as a comma-separated list. '
       'Default: "%(default)s".')
args = parser.parse_args()

# Read arguments.
try:
  allowed_flags = map(int, args.flags.split(','))
except ValueError:
  sys.stderr.write('Error: Invalid flags "{}". Must be a comma-separated list of integers.\n'
                   .format(args.flags))
  sys.exit(1)

if args.output:
  fasta_path = args.output
else:
  fasta_path = args.bam+'.major.fa'
samfile = pysam.Samfile(args.bam, 'rb')

def select_read(read):
  """Return the read if it passes the thresholds, None otherwise."""
  if (read.alignment.mapq >= args.min_mapq and
      read.alignment.flag in allowed_flags):
    return read
  else:
    return False

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

sequence = ['N'] * args.length

# For each position in the reference:
for pileupcolumn in samfile.pileup(args.chrom, 0, args.length, stepper='all', max_depth=10000000,
                                   mask=False):
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
      if ord(alignment.qual[pos_in_read])-33 >= args.min_phred:
        # Then count it in the list of bases at this position
        bases_in_position.append(alignment.seq[pos_in_read])
  # Find the major allele and add it to the sequence
  if bases_in_position:
    sequence[pos_in_ref] = get_major(bases_in_position)
  # If there are no reads passing thresholds, the base will be "N".

# Write consensus sequence to output FASTA file
with open(fasta_path, 'w+') as fastafile:
  printseq(sequence, fastafile, args.chrom)
