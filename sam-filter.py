#!/usr/bin/env python
# Notes:
# Test data has shown an example of a deletion occurring beyond the border of
# the chromosome. Just check that it's in a valid coordinate range.
# - occurred in M01368:8:000000000-A3GHV:1:1108:18888:4511 (sample R19S5)
from __future__ import division
import re
import os
import sys
from optparse import OptionParser

USAGE = "USAGE: %prog [options] file.sam"
DESCRIPTION = """Slightly more advanced SAM filtering than SAMTools. E.g. CIGAR
parsing, NM tag, etc.
Give "-" as the input filename to read from stdin."""
EPILOG = """"""
OPT_DEFAULTS = {'region':'', 'event_types':'', 'clipping':0.0, 'coord':0,
  'sam_output':False, 'no_header':False, 'test_output':False}

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
for opt in sys.argv[1:]:
  if opt == '--debug' or re.search(r'^-\w*d', opt):
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
    help='Filter for CIGAR string events. Reads will only pass this filter if '
      +'they have one of the given CIGAR events that *starts* at the '
      +'coordinate given with -c (if no coordinate is given, any reads with '
      +'the event anywhere will pass. Specify the event types by their CIGAR '
      +'characters, e.g. "ID" for insertions and deletions (currently the only '
      +'supported ones).')
  parser.add_option('-C', '--clipping', dest='clipping', type='float',
    default=OPT_DEFAULTS.get('clipping'),
    help='Filter out reads which have more than this much clipping, as a '
      +'percentage of their length.')
  parser.add_option('-I', '--indel-output', dest='indel_output',
    action='store_true', default=OPT_DEFAULTS.get('indel_output'),
    help='Just print information on found indels passing the filters, instead '
    +'of SAM output. This is a legacy output. If it is enabled, all filters '
    +'except the CIGAR event types one will be ignored.')
  parser.add_option('-H', '--no-header', dest='no_header',
    action='store_true', default=OPT_DEFAULTS.get('no_header'),
    help='Don\'t include the SAM header.')
  parser.add_option('-T', '--test-output', dest='test_output',
    action='store_true', default=OPT_DEFAULTS.get('test_output'),
    help='Generate original test output format on the data.')
  parser.add_option('-d', '--debug', dest='debug', action='store_true',
    help='Turn on debug mode.')

  (options, arguments) = parser.parse_args()

  indel_output = options.indel_output
  test_output = options.test_output
  event_types = options.event_types
  clipping = options.clipping
  no_header = options.no_header
  region = options.region
  coord = options.coord

  if not arguments:
    parser.print_help()
    fail("\nProvide a filename (or \"-\") as an argument")
  else:
    infilename = arguments[0]

  if infilename == '-':
    infilehandle = sys.stdin
  else:
    infilehandle = open(infilename, 'r')

  for line in infilehandle:
    passed = True

    line_stripped = line.strip()
    # empty lines
    if not line:
      continue
    # header lines
    if line_stripped[0] == '@' and line_stripped[3] == "\t":
      if not (no_header or indel_output):
        sys.stdout.write(line)
      continue

    fields = line_stripped.split("\t")
    if len(fields) < 11:
      fail("Format Error: Too few columns in SAM file. Line:\n"+line)

    name  = fields[0]
    cigar = fields[5]
    seq   = fields[9]
    pos   = int(fields[3])
    if test_output:
      print "\n"+name+"\n"+str(pos)+": "+cigar

    # CIGAR events filter
    if event_types:
      passed_filter = False
      indels = parse_cigar(cigar, pos, seq=seq, test_output=test_output)
      for indel in indels:
        if indel[0] in event_types and (not coord or coord == indel[1]):
          passed_filter = True
          if indel_output:
            print indel
          else:
            break
      # if indel_output, skip rest of main loop to avoid printing SAM output
      if passed_filter and indel_output:
        continue
      passed = passed and passed_filter

    # clipping filter
    if clipping:
      passed_filter = clip_pass(cigar, clipping)
      passed = passed and passed_filter

    if passed:
      sys.stdout.write(line)

  if isinstance(infilehandle, file):
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


def clip_pass(cigar, threshold):
  """Find whether the read passes the soft clipping threshold.
  threshold is a maximum soft-clipping percentage (float)."""

  clip_length = 0
  total_length = 0
  for move in re.findall(r'\d*[MIDNSHPX=]', cigar):
    op = move[-1]
    if len(move) > 1:
      op_length = int(move[:-1])
    else:
      op_length = 1
      sys.stderr.write("Warning: CIGAR op with no length (assumed 1): "+cigar)
    if op == 'M':
      total_length += op_length
    elif op == 'I':
      total_length += op_length
    elif op == 'D':
      pass
    elif op == 'N':
      pass
    elif op == 'S':
      total_length += op_length
      clip_length += op_length
    elif op == 'H':
      total_length += op_length
      clip_length += op_length
    elif op == 'P':
      pass

  # print "clipped:%6.2f (%s)" % (100*clip_length/total_length,
  #   ' '.join(re.findall(r'\d*[MIDNSHPX=]', cigar)))

  if 100*clip_length/total_length > threshold:
    return False
  else:
    return True


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()