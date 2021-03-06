#!/usr/bin/env python
#TODO: Warn when site spans the edge of an alignment block
#TODO: Reverse complement sequence
from __future__ import division
import sys
import argparse
import lavreader
import lavintervals

EXPECTED_VERSIONS = {'lavreader':'0.71', 'lavintervals':'0.81'}

OPT_DEFAULTS = {'tsv':True, 'slop':20}
USAGE = "%(prog)s [options] align.lav (sites.tsv|-s chr:coord)"
DESCRIPTION = """Default input format: the tsv output of inspect-reads.py.
Prints new version to stdout. If a site lies in an unknown region (does not
align to the reference), it will be left unaltered (unless -N is specified).
If a site lies in multiple contigs, the alignment data from the largest one
will be used to compute the new coordinates."""
EPILOG = """"""

def main():
  version_check(EXPECTED_VERSIONS)

  parser = argparse.ArgumentParser(description=DESCRIPTION)#, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('lavpath', metavar='align.lav',
    help='The LASTZ alignment, in LAV format. Should be produced with the '
      'command "$ lastz reference.fa assembly.fa". This was developed for '
      'LASTZ version 1.02.00.')
  parser.add_argument('sites_file', metavar='sites.tsv', nargs='?',
    help='A file containing the sites to convert. By default, should be in '
      'a tsv returned by inspect-reads.py. Can be omitted, if sites are given '
      'on command line.')
  parser.add_argument('-t', '--tsv', action='store_true',
    help='Input sites file is in the tsv format produced by inspect-reads.py '
      '(default mode).')
  parser.add_argument('-v', '--vcf', action='store_true',
    help='Input sites file is a VCF.')
  parser.add_argument('--string', action='store_true',
    help='Output the lifted-over sites in UCSC format (see --sites for format '
      'details).')
  parser.add_argument('-s', '--sites', nargs='*', metavar='SITE',
    help='The sites to convert, given as a series of UCSC-format coordinates. '
      'Allows giving the sites directly as command-line arguments, instead of '
      'in a file. Example: "-s chr3:513 chrM:2456 21:4090".')
  parser.add_argument('-S', '--slop',
    help='The "slop" value used for asm-curator.py. Default: %(default)s.')
  parser.add_argument('-N', '--null-note', metavar='NOTE',
    help='If a site lies in an unknown region (does not align to the '
      'reference), append this string to the chromosome name.')

  args = parser.parse_args()
  if not (args.sites or args.sites_file):
    parser.print_help()
    fail('Error: Must provide sites in either a file or command-line arguments.')

  sites = read_sites(args)
  lav = lavreader.LavReader(args.lavpath)

  # Only use the contigs kept by asm-curator.py, for simplicity.
  #TODO: allow skipping curation
  contigs = get_kept_contigs(lav, slop=args.slop)
  # Get the coordinate conversion coefficients for every gap-free block.
  table = lavintervals.hits_to_conv_table(lav, contigs=contigs, tuples=True)

  # For each site, compare it to every block to check if it's contained in it.
  # Use the longest one.
  for (s_chrom, s_begin, s_end, string) in sites:
    if s_begin is None:
      # if no start coordinate, just print the input line back out
      sys.stdout.write(string)
      continue
    # Get a list of all blocks that contain the site.
    containing_blocks = []
    for block in table:
      (b_chrom1, b_begin, b_end) = block[:3]
      if s_chrom == b_chrom1 and b_begin <= s_begin <= b_end:
        containing_blocks.append(block)
    # If there are no blocks that contain the site, there is no reference
    # sequence corresponding to this region. Just print the old site unedited.
    if not containing_blocks:
      if args.null_note:
        raise NotImplementedError('--null-note not yet implemented.')
      else:
        sys.stdout.write(string)
      continue
    # Find the longest block among the hits.
    longest_begin = 0
    longest_end = 0
    for block in containing_blocks:
      (b_begin, b_end) = block[1:3]
      if b_end - b_begin > longest_end - longest_begin:
        longest_block = block
        (longest_begin, longest_end) = longest_block[1:3]
    (b_chrom1, b_begin, b_end, b_chrom2, b_strand, b_offset, b_id, b_score) = longest_block
    # Do the actual conversion.
    # get the chromosome name
    s_chrom = b_chrom2
    # calculate the new coordinate, using the conversion coefficients
    s_begin = s_begin * b_strand + b_offset
    if s_end:
      s_end = s_end * b_strand + b_offset
    site = (s_chrom, s_begin, s_end, string)
    sys.stdout.write(edit_line(site, args, strand=b_strand))


def get_kept_contigs(lav, slop=20):
  """Return the contigs which would be retained by asm-curator.py.
  Uses the same procedure as in that script."""
  contigs = set()
  #TODO: replace with reading retained contigs from FASTA
  intervals = lavintervals.alignments_to_intervals(lav)
  all_overlaps = lavintervals.get_all_overlaps(intervals)
  all_overlaps = lavintervals.discard_redundant(all_overlaps, slop=slop)
  # remove discarded intervals
  for interval in intervals.keys():
    if interval in all_overlaps:
      hit = intervals[interval].parent
      contig = hit.query['name'].split()[0]
      contigs.add(contig)
  return contigs


#TODO: Join into one function (they're mostly the same, it ends up)
#TODO: Compute site end
def read_sites(args):
  """Read in a sites file (or string), return a generator.
  "args" must have attributes 'sites', 'tsv', and 'vcf' to input source, and
  'sites_file' if source is a file.
  For each line, the generator returns a list of 4 values:
  0 (chrom): the chromosome of the variant (str)
  1 (begin): the start coordinate of the variant (int)
  2 (end):   not yet implemented (always None)
  3 (line):  the raw input line (str)
  For non-data lines (headers, etc), the first 3 values will be None.
  """
  if args.sites:
    return parse_sites_string(args.sites)
  elif args.vcf:
    return read_vcffile(args.sites_file)
  elif args.tsv:
    return read_tsvfile(args.sites_file)
  else:
    fail("Error: Unable to determine input filetype.")


def parse_sites_string(sites_strs):
  for site_str in sites_strs:
    try:
      (chrom, coord) = site_str.split(':')
      yield [chrom, int(coord), None, site_str+'\n']
    except ValueError:
      fail('Error: Invalid site string "'+site_str+'"')


def read_tsvfile(filepath):
  """Read in a tsv file produced from inspect-reads.py.
  Return data structure described in read_sites_file().
  """
  with open(filepath, 'r') as filehandle:
    for raw_line in filehandle:
      line = raw_line.strip()
      # don't process empty or commented lines (includes header line)
      if not line or line.startswith('#'):
        yield [None, None, None, raw_line]
        continue
      fields = raw_line.split('\t')
      try:
        yield [fields[1], int(fields[2]), None, raw_line]
      except (IndexError, ValueError):
        fail("Error: Input file "+filepath+" not in the tsv format expected "
          "from inspect-reads.py.")


def read_vcffile(filepath):
  with open(filepath, 'r') as filehandle:
    for raw_line in filehandle:
      line = raw_line.strip()
      if not line or line.startswith('#'):
        yield [None, None, None, raw_line]
        continue
      fields = raw_line.split('\t')
      try:
        yield [fields[0], int(fields[1]), None, raw_line]
      except (IndexError, ValueError):
        fail("Error: Input file "+filepath+" not in VCF.")


#TODO: Join into one function (they're mostly the same, it ends up)
def edit_line(site, args, strand=1):
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
  if args.string or args.sites:
    return edit_site_string(site, strand=strand)
  elif args.vcf:
    return edit_vcfline(site, strand=strand)
  elif args.tsv:
    return edit_tsvline(site, strand=strand)
  else:
    fail("Error: Unrecognized filetype.")


def edit_tsvline(site, strand=1):
  """See edit_line() for interface."""
  chrom = site[0]
  begin = site[1]
  # end = site[2]
  raw_line = site[3]
  # do not alter empty/header lines
  if begin is None:
    return raw_line
  fields = raw_line.split('\t')
  if len(fields) < 3:
    fail('Error: Raw line not in the tsv format expected from '
      'inspect-reads.py: "'+raw_line+'"')
  # adjust coordinates for reverse complement indels
  if strand == -1 and len(fields) > 3 and fields[3] in 'ID':
    begin -= 2
  if chrom is not None:
    fields[1] = chrom
  # preserve newline if it's the last field
  if fields[2][-2:] == '\r\n':
    begin = str(begin)+fields[1][-2:]
  elif fields[2][-1] in '\r\n':
    begin = str(begin)+fields[1][-1]
  fields[2] = begin
  #TODO: make sure newline is preserved
  return '\t'.join(map(str, fields))


def edit_vcfline(site, strand=1):
  """See edit_line() for interface."""
  chrom = site[0]
  begin = site[1]
  end = site[2]
  raw_line = site[3]
  # do not alter empty/header lines
  if begin is None:
    return raw_line
  fields = raw_line.split('\t')
  try:
    if chrom is not None:
      fields[0] = chrom
    if begin is not None:
      fields[1] = str(begin)
    if strand == -1:
      fields[3] = 'r'+fields[3]
      fields[4] = 'r'+fields[4]
  except (IndexError, TypeError):
    fail('Error: Raw line not in VCF: "'+raw_line+'"')
  #TODO: make sure newline is preserved
  return '\t'.join(map(str, fields))


def edit_site_string(site, strand=1):
  chrom = site[0]
  begin = site[1]
  end = site[2]
  site_str = site[3]
  # do not alter empty/header lines
  if begin is None:
    return site_str
  fields = site_str.split(':')
  try:
    if chrom is not None:
      fields[0] = chrom
    if begin is not None:
      fields[1] = str(begin)
    if end is not None:
      fields[2] = str(end)
  except (IndexError, TypeError):
    fail('Error: Invalid site string: "'+site_str+'"')
  #TODO: make sure newline is preserved
  if end is None:
    return chrom+':'+str(begin)+'\n'
  else:
    return chrom+':'+str(begin)+'-'+str(end)+'\n'


def version_check(expected):
  actual = {}
  for module_name in expected:
    module = sys.modules[module_name]
    for version_name in ['version', 'VERSION', '__version__']:
      if version_name in dir(module):
        actual[module_name] = getattr(module, version_name)
  for module_name in actual:
    assert actual[module_name] == expected[module_name], (
      "Wrong version of "+module_name+". Expected: "+expected[module_name]
      +", actual: "+actual[module_name]
    )


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
