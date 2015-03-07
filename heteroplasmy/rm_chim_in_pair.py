#!/usr/bin/env python
import os
import pysam
import argparse

OPT_DEFAULTS = {'length':100, 'region':'chrM', 'bounds':(600,16000)}
DESCRIPTION = """Remove chimeric reads, and reads under a minimum length from an
an alignment. Preserves read pairs by removing either both or neither."""
USAGE = """%(prog)s [-r region] input.bam [output.bam]"""

parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.set_defaults(**OPT_DEFAULTS)
parser.add_argument('input', metavar='input.bam',
  help='The input alignment.')
parser.add_argument('output', metavar='output.bam', nargs='?',
  help='The output alignment (optional). If not given, the output filename '
    'will be the input filename prepended with "dechim.rlen.".')
parser.add_argument('-l', '--length', metavar='length', type=int,
  help='Minimum read length. Default: %(default)s.')
parser.add_argument('-r', '--region', metavar='chrom',
  help='The name of the target sequence (chromosome). Default: %(default)s.')
parser.add_argument('-b', '--bounds', metavar='coord', type=int, nargs=2,
  help='The bounds of the chimera-free region of the target. Reads whose '
    'chimeric partner is located outside this region will still pass the '
    'filter. Intended for circular genomes and plasmids where reads will map '
    'partially to the start and end at the same time. Default: %(default)s.')
parser.add_argument('-m', '--margin', type=int,
  help='The size of the chimera-free regions at each end of the target. Using '
    'this and --chr-length is an alternative to using --bounds.')
parser.add_argument('-L', '--chr-length', type=int,
  help='The length of the target sequence.')
args = parser.parse_args()

bounds = args.bounds
if args.margin and args.chr_length:
  bounds = (args.margin, args.chr_length-args.margin)
elif args.margin or args.chr_length:
  raise Exception('--margin and --chr-length must be used together.')

if not args.output:
  dirpath, filename = os.path.split(args.input)
  args.output = os.path.join(dirpath, 'dechim.rlen.'+filename)


def check_chim(read):
  tags = dict(read.tags)
  if 'SA' in tags.keys():
    sa = tags['SA'].split(';')[:-1]
    if len(sa) == 1:
      (chrom, pos, strand, cigar, mapq, nm) = sa[0].split(',')
      if (chrom == args.region and
          (int(pos) <= bounds[0] or int(pos) >= bounds[1])):
        return read
      else:
        pass
    else:
      pass
  else:
    return read

sam = pysam.Samfile(args.input, 'rb')
out = pysam.Samfile(args.output, 'wb', template=sam)
for read in sam:
  try:
    read1 = read
    read2 = sam.next()
    if (check_chim(read1) and check_chim(read2) and
        read1.rlen >= args.length and read2.rlen >= args.length):
      out.write(read1)
      out.write(read2)
    else:
      pass
  except StopIteration:
    sam.close()
    out.close()

