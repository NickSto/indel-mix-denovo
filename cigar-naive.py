#!/usr/bin/env python
# Notes:
# Test data has shown an example of a deletion occurring beyond the border of
# the chromosome. Just check that it's in a valid coordinate range.
# - occurred in M01368:8:000000000-A3GHV:1:1108:18888:4511 (sample R19S5)
import re
import os
import sys
from optparse import OptionParser

USAGE = "USAGE: %prog [options] file.sam"
DESCRIPTION = """A quick-and-dirty script to parse the CIGAR strings in a SAM
file.
Give "-" as the input filename to read from stdin."""
EPILOG = """"""
OPT_DEFAULTS = {'region':'', 'event_types':'ID', 'coord':0, 'sam_output':False,
  'no_header':False, 'test_output':False}

TESTS = [
  '251M',       # M
  '3S248M',     # SM
  '4M2I245M',   # MIM
  '166M1D85M',  # MDM
  '179M62I2M2D2M6S',      # MIMDMS
  '199M1D2M2I8M2D2M38S',  # MDMIMDMS
]

# try to set debug as a global
debug = False
if '-d' in sys.argv[1:] or '--debug' in sys.argv[1:]:
  debug = True
for opt in sys.argv[1:]:
  if re.match(r'^-[a-zA-Z]*d', opt):
    debug = True

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-r', '--region', dest='region',
    default=OPT_DEFAULTS.get('region'),
    help='***NOT IMPLEMENTED YET*** The region to restrict the search to (UCSC '
      +'coordinate format)')
  parser.add_option('-c', '--coord', dest='coord', type='int',
    default=OPT_DEFAULTS.get('coord'),
    help='Give a specific coordinate to restrict the indel search to.')
  parser.add_option('-e', '--event-types', dest='event_types',
    default=OPT_DEFAULTS.get('event_types'),
    help=('The CIGAR operations to output: give a string of letters for each '
      +'you want to include, e.g. "ID" for insertions and deletions (currently '
      +'the only supported ones). Defaults: "%default"'))
  parser.add_option('-s', '--sam-output', dest='sam_output',
    action='store_const', const=not(OPT_DEFAULTS.get('sam_output')),
    default=OPT_DEFAULTS.get('sam_output'),
    help='When a read with a matching indel is found, print the original line '
      +'from the SAM file instead of the indel data.')
  parser.add_option('-H', '--no-header', dest='no_header',
    action='store_const', const=not(OPT_DEFAULTS.get('no_header')),
    default=OPT_DEFAULTS.get('no_header'),
    help='If using SAM output (-s), don\'t include the header.')
  parser.add_option('-T', '--test-output', dest='test_output',
    action='store_const', const=not(OPT_DEFAULTS.get('test_output')),
    default=OPT_DEFAULTS.get('test_output'),
    help='Generate original test output format on the data.')
  parser.add_option('-d', '--debug', dest='debug', action='store_true',
    help=('Turn on debug mode.'))

  (options, arguments) = parser.parse_args()

  test_output = options.test_output
  event_types = options.event_types
  sam_output = options.sam_output
  no_header = options.no_header
  region = options.region
  coord = options.coord

  if not arguments:
    parser.print_help()
    fail("Provide a filename (or \"-\") as an argument")
  else:
    infilename = arguments[0]

  if infilename == '-':
    infilehandle = sys.stdin
  else:
    infilehandle = open(infilename, 'r')

  for line_orig in infilehandle:

    line = line_orig.strip()
    if not line:
      continue
    if line[0] == '@' and line[3] == "\t":
      if sam_output and not no_header:
        sys.stdout.write(line_orig)
      continue

    fields = line.split("\t")

    name  = fields[0]
    cigar = fields[5]
    seq   = fields[9]
    pos   = int(fields[3])
    if test_output:
      print "\n"+name+"\n"+str(pos)+": "+cigar
    indels = parse_cigar(cigar, pos, seq=seq, test_output=test_output)
    for indel in indels:
      if indel[0] in event_types:
        if not coord or coord == indel[1]:
          if sam_output:
            sys.stdout.write(line_orig)
            break
          else:
            print indel

  if infilehandle is not sys.stdin:
    infilehandle.close()


def parse_cigar(cigar, pos, seq='', test_output=False):
  """Returns a list of indels found in the alignment. Each indel is a tuple:
  (type, reference_coordinate, length)
  type: 'I' for insertion or 'D' for deletion
  reference_coordinate: the 1-based coordinate of the base preceding the indel
  length: its length"""
  indels = []

  ref_pos = pos
  alt_pos = 0
  for move in re.findall(r'[^A-Z]*[A-Z]', cigar):
    match = re.search(r'^(\d+)([MIDSHPN])$', move)
    if match:
      length = int(match.group(1))
      op = match.group(2)
    else:
      match = re.search(r'^([MIDSHPN])$', move)
      if match:
        length = 1
        op = match.group(1)
      else:
        fail("Error: CIGAR contains a movement in unrecognized format.\n"
          +"CIGAR: "+cigar+"\tmovement: "+move)
    if op == 'M':
      ref_pos += length
      alt_pos += length
    elif op == 'I':
      # ref_pos-1 is the (1-based) base to the left of the insertion
      # if ref_pos-1 > 0 and ref_pos < ref_end:
      if test_output:
        print ("insertion after  "+str(ref_pos-1)+": "+seq[alt_pos-4:alt_pos]+"{"
          +seq[alt_pos:alt_pos+length]+"}"+seq[alt_pos+length:alt_pos+length+4])
      indels.append(('I', ref_pos-1, length))
      alt_pos += length
    elif op == 'D':
      # if ref_pos-1 > 0 and ref_pos+length < ref_end:
      if test_output:
        print ("deletion between "+str(ref_pos-1)+" and "+str(ref_pos+length)+": "
          +seq[alt_pos-4:alt_pos]+"["+("-"*length)+"]"
          +seq[alt_pos:alt_pos+4])
      indels.append(('D', ref_pos-1, length))
      ref_pos += length
    elif op == 'N':
      if test_output:
        print "N from "+str(ref_pos)+" to "+str(ref_pos+length)
      ref_pos += length
    elif op == 'S':
      alt_pos += length
    elif op == 'H':
      pass
    elif op == 'P':
      pass

  return indels


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()