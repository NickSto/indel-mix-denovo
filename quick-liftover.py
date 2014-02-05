#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import lavreader
import lavintervals

OPT_DEFAULTS = {'tsv':True}
USAGE = "USAGE: %(prog)s [options] align.lav sites.tsv"
DESCRIPTION = """Default input format: the tsv output of inspect-reads.py."""
EPILOG = """"""

SLOP = 20

def main():

  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('lavpath', metavar='align.lav',
    help="""""")
  parser.add_argument('sitespath', metavar='sites.tsv',
    help="""""")
  parser.add_argument('-t', '--tsv', action='store_true',
    help="""Input data (second file) is in the tsv format produced by
      inspect-reads.py (default).""")
  parser.add_argument('-v', '--vcf', action='store_true',
    help="""Input data (second file) is a VCF.""")
  parser.add_argument('-f', '--float', type=float,
    help="""default: %(default)s""")

  args = parser.parse_args()

  sitesreader = read_sitesfile(sitespath, args)
  lav = lavreader.LavReader(args.lavpath)

  contigs = get_kept_contigs(lav)
  table = blocks_to_conv_table(lav, contigs=contigs)

  # for each site:
  # just compare to each interval to check if it's contained
  # use the longest one
  for (site_chrom, site_begin, site_end, line) in sitesreader:
    if begin is None:
      print line
    for (block_chrom, block_begin, block_end, refname, strand, offset) in table:
      if (site_chrom == block_begin and
          block_begin <= site_begin <= block_end):
        site_chrom = refname
        site_begin = site_begin * strand + offset
        site_end = site_end * strand + offset
        
    print edit_line(site, args)


def get_kept_contigs(lav):
  contigs = set()
  #TODO: replace with reading retained contigs from FASTA
  intervals = lavintervals.alignments_to_intervals(lav)
  all_overlaps = lavintervals.get_all_overlaps(intervals)
  all_overlaps = lavintervals.discard_redundant(all_overlaps, slop=SLOP)
  # remove discarded intervals
  for interval in intervals.keys():
    if interval in all_overlaps:
      hit = intervals[interval].parent
      contigs.add(hit.query['name'])
  return contigs


def read_sitesfile(filepath, args):
  """Read in a sites file, return a generator.
  "args" must have bool attributes args.tsv and args.vcf to determine filetype.
  For each line, the generator returns a tuple of 4 values:
  0 (chrom): the chromosome of the variant (str)
  1 (begin): the start coordinate of the variant (int)
  2 (end):   not yet implemented (always None)
  3 (line):  the raw input line (str)
  For non-data lines (headers, etc), the first 3 values will be None.
  """
  if args.tsv:
    return read_tsvfile(filepath)
  elif args.vcf:
    fail("Error: VCF input not yet supported.")
  else:
    fail("Error: Must specify the sites input filetype.")


def read_tsvfile(filepath):
  """Read in a tsv file produced from inspect-reads.py.
  Return data structure described in read_sitesfile().
  """
  with open(filepath, 'r') as filehandle:
    for raw_line in filehandle:
      line = raw_line.strip()
      # don't process empty or commented lines (includes header line)
      if not line or line.startswith('#'):
        yield (None, None, None, raw_line)
      fields = line.split('\t')
      try:
        yield (fields[0], int(fields[1]), None, raw_line)
      except (IndexError, ValueError):
        fail("Error: Input file "+filepath+" not in the tsv format expected "
          +"from inspect-reads.py.")


def edit_line(site, args):
  """Edit a raw line from an input file to reflect the converted coordinates.
  I.e. It replaces the chromosome and coordinate data in the raw line (the last
  value in "site") with the first three values in "site". Any None value will
  leave the original field unedited.
  "site" should be a site/line in the same format produced from the appropriate
    file reading function (e.g. read_tsvstats()).
  "args" should have the following bool attributes indicating filetype:
    'tsv' - tsv from inspect-reads.py
    'vcf' - VCF (not yet implemented)
  Output is the final line, ready for printing (includes newline).
  """
  if args.tsv:
    return edit_tsvline(site)
  elif args.vcf:
    fail("Error: VCF input not yet supported.")
  else:
    fail("Error: Unrecognized filetype.")


def edit_tsvline(site):
  """See edit_line() for interface."""
  chrom = site[0]
  begin = site[1]
  end = site[2]
  raw_line = site[3]
  # do not alter empty/header lines
  if begin is None:
    return raw_line
  fields = raw_line.strip().split('\t')
  try:
    if chrom is not None:
      fields[0] = chrom
    if begin is not None:
      fields[1] = str(begin)
  except IndexError:
    fail("Error: Raw line not in the tsv format expected from inspect-reads.py")
  #TODO: make sure newline is preserved
  return '\t'.join(fields)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
