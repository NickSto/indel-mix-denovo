#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import argparse

LABELS = ('Family', 'Sample', 'Chrom', 'Coord', 'Type', 'Alt', 'Covg', 'Reads', 'Freq', 'Unmap',
          'Improp', 'Fwd', 'First', 'Sndary', 'Dup', 'Mapq0', 'Mapq20', 'Mapq30', 'MapqMax',
          'MapqMaxReads', 'SBias', 'MBias', 'Flags', 'PosDist', 'Seq', 'Context')

ARG_DEFAULTS = {'input':sys.stdin, 'covg':1000, 'freq':1.0, 'strand':1.0, 'mate':1.0,
                'excluded':'chrM:302-310,chrM:16183-16192', 'pos_method':'range', 'pos_thres':3}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='vars_asm.tsv', type=argparse.FileType('r'), nargs='?',
    help='Output of inspect-reads.py')
  parser.add_argument('-c', '--covg', type=int,
    help='Coverage threshold. Variants with less than this depth of coverage will be filtered out. '
         'Default: %(default)s')
  parser.add_argument('-f', '--freq', type=float,
    help='Minor allele frequency threshold, in percent. Variants with a lower MAF will be filtered '
         'out. Default: %(default)s%%')
  parser.add_argument('-s', '--strand', type=float,
    help='Strand bias threshold. Variants with a higher strand bias will be filtered out. '
         'Default: %(default)s')
  parser.add_argument('-m', '--mate', type=float,
    help='Mate bias threshold. Variants with a higher mate bias will be filtered out. '
         'Default: %(default)s')
  parser.add_argument('-e', '--excluded',
    help='Regions to exclude. Variants in these regions will be filtered out. '
         'Format for excluded: comma-separated regions. Format for regions: '
         '"[chrom]:[start]-[end]", e.g. "302-310", "chrM:302-310", "SC8-ch:302-310", '
         '"1:1507-1510", "1-12:1002-1100". [chrom] is optional. If omitted, the region will be '
         'assumed to apply to every chromosome. Default: %(default)s')
  parser.add_argument('-M', '--pos-method', choices=('range', 'r', 'rsquared'))
  parser.add_argument('-t', '--pos-thres', type=float,
    help='Default: %(default)s')
  parser.add_argument('-r', '--max-r', type=float,
    help='Default: %(default)s')
  parser.add_argument('-D', '--debug', action='store_true')

  args = parser.parse_args(argv[1:])

  excluded_regions = parse_excluded(args.excluded)

  labels = None
  header = True
  for line in args.input:
    fields = line.rstrip('\r\n').split('\t')
    # At the header (first line), parse the labels list.
    if header:
      assert is_header(fields), 'First line should be a header with field labels.'
      labels = [label.lower() for label in fields]
      labels[0] = labels[0].lstrip('#')
      sys.stdout.write(line)
      header = False
    # Skip header lines interspersed in-between data lines.
    if is_header(fields):
      continue
    # Make a dict with the variant data. Try to parse each value as an int, float, etc.
    variant = {}
    for label, field in zip(labels, fields):
      if label == 'posdist':
        if field == '.' or field == '':
          value = None
        else:
          value = [int(count) for count in field.split(',')]
      else:
        try:
          value = int(field)
        except ValueError:
          try:
            value = float(field)
          except ValueError:
            value = field
      variant[label] = value
    # Apply filters. Start with faster filters first.
    passed = True
    if variant['freq'] < args.freq:
      passed = False
    elif variant['covg'] < args.covg:
      passed = False
    elif variant['sbias'] > args.strand:
      passed = False
    elif variant['mbias'] > args.mate:
      passed = False
    elif is_excluded(variant, excluded_regions):
      passed = False
    elif is_position_biased(variant, args.pos_method, args.pos_thres):
      passed = False
    # Print the variant, if it passed.
    if passed:
      sys.stdout.write(line)
    elif args.debug:
      pos_dist = variant['posdist']
      bins = len(pos_dist)
      expected = variant['reads']/bins
      pos_dist_ratio = ' '.join(['{:4.1f}'.format(n/expected) for n in pos_dist])
      pos_dist_raw   = ' '.join(['{:3d}'.format(n) for n in pos_dist])
      print('{coord:9s}{type:5s}{alt:10s}{covg:9s}{freq:9s}{sbias:10s}{mbias:10s}{reads:5s}|  '
            '{raw:43s} {ratio:43s}'.format(ratio=pos_dist_ratio, raw=pos_dist_raw, **variant))


def is_header(fields):
  """Try to heuristically determine whether the line is a header.
  Returns True if we can say confidently it's a header line, False if we can say confidently it's
  not, and None if we can't determine. Since we'll use this on every line of thousands, the tests
  are organized to be as quick as possible in the typical case."""
  if len(fields) < 11:
    return None
  # Test a column which should usually be a number in data lines and never a number in header lines.
  try:
    float(fields[8])
    return False
  except ValueError:
    pass
  first_field = fields[0]
  # An explicitly commented line is a header.
  if first_field.startswith('#'):
    return True
  # The first field in a header is usually these two (and never these in data lines).
  if first_field.lower() == 'sample' or first_field.lower() == 'family':
    return True
  # Fallback 1: There should never be a number in a header line. If we find one, it's a data line.
  for field in fields:
    try:
      float(field)
      return False
    except ValueError:
      pass
  # Fallback 2: Just test whether any of the known labels is in the line.
  for label in LABELS:
    if label in fields:
      return True
  for label in LABELS:
    if label.lower() in fields:
      return True


def parse_excluded(excluded_str):
  # Regions: "302-310", "chrM:302-310", "SC8-ch:302-310", "1:1507-1510", "1-12:1002-1100"
  regions = []
  region_strs = excluded_str.split(',')
  for region_str in region_strs:
    region = {'chrom':None}
    fields = region_str.split(':')
    if len(fields) == 1:
      region['chrom'] = None
      range_str = region_str
    else:
      assert len(fields) == 2, 'Excluded region "'+region_str+'" contains too many colons.'
      region['chrom'] = fields[0]
      range_str = fields[1]
    fields = range_str.split('-')
    assert len(fields) == 2, 'Excluded coordinate range "'+range_str+'" must contain one dash.'
    try:
      region['start'] = int(fields[0])
    except ValueError:
      sys.stderr.write('Excluded coordinate "'+fields[0]+'" not a number.')
      raise
    try:
      region['end'] = int(fields[1])
    except ValueError:
      sys.stderr.write('Excluded coordinate "'+fields[1]+'" not a number.')
      raise
    regions.append(region)
  return regions


def is_excluded(variant, excluded_regions):
  for region in excluded_regions:
    if ((region['chrom'] is None or region['chrom'] == variant['chrom']) and
        region['start'] <= variant['coord'] <= region['end']):
      return True
  return False


def is_position_biased(variant, method, thres):
  # Filter by distribution of variant position in the read.
  pos_dist = variant['posdist']
  total = variant['reads']
  bins = len(pos_dist)
  expected = total/bins
  min_count = expected/thres
  max_count = thres*expected
  r = 0
  for count in pos_dist:
    if method == 'range':
      if count < min_count or count > max_count:
        return True
    elif method == 'r':
      r += abs(count - expected)
  if method == 'r':
    if r/total > thres:
      return True


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
