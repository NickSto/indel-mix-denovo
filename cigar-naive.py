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
OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'test_output':False}

TESTS = [
  '251M',       # M
  '3S248M',     # SM
  '4M2I245M',   # MIM
  '166M1D85M',  # MDM
  '179M62I2M2D2M6S',      # MIMDMS
  '199M1D2M2I8M2D2M38S',  # MDMIMDMS
]

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-s', '--str', dest='str',
    default=OPT_DEFAULTS.get('str'), help='default: %default')
  parser.add_option('-i', '--int', dest='int', type='int',
    default=OPT_DEFAULTS.get('int'), help='')
  parser.add_option('-f', '--float', dest='float', type='float',
    default=OPT_DEFAULTS.get('float'), help='')
  parser.add_option('-T', '--test-output', dest='test_output',
    action='store_const', const=not OPT_DEFAULTS.get('test_output'),
    default=OPT_DEFAULTS.get('test_output'),
    help='Generate original test output format on the data.')

  (options, arguments) = parser.parse_args()

  test_output = options.test_output

  if not arguments:
    parser.print_help()
    fail("Provide a filename (or \"-\") as an argument")
  else:
    infilename = arguments[0]

  if infilename == '-':
    infilehandle = sys.stdin
  else:
    infilehandle = open(infilename, 'r')

  for line in infilehandle:

    line = line.strip()
    if not line:
      continue
    if line[0] == '@' and line[3] == "\t":
      continue # header line

    fields = line.split("\t")

    name  = fields[0]
    cigar = fields[5]
    seq   = fields[9]
    pos   = int(fields[3])
    if test_output:
      print "\n"+name+"\n"+str(pos)+": "+cigar
    parse_cigar(cigar, seq, pos, test_output)

  if infilehandle is not sys.stdin:
    infilehandle.close()


def parse_cigar(cigar, seq, pos, test_output):

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
      alt_pos += length
    elif op == 'D':
      # if ref_pos-1 > 0 and ref_pos+length < ref_end:
      if test_output:
        print ("deletion between "+str(ref_pos-1)+" and "+str(ref_pos+length)+": "
          +seq[alt_pos-4:alt_pos]+"["+("-"*length)+"]"
          +seq[alt_pos:alt_pos+4])
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


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()