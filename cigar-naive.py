#!/usr/bin/env python
# Notes:
# Test data has shown an example of a deletion occurring beyond the border of
# the chromosome. Just check that it's in a valid coordinate range.
# - occurred in M01368:8:000000000-A3GHV:1:1108:18888:4511 (sample R19S5)
import re
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'test':False}
USAGE = "USAGE: %prog [options] file.sam"
DESCRIPTION = """Give "-" as the input filename to read from stdin."""
EPILOG = """"""

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
  parser.add_option('-t', '--test', dest='test', action='store_const',
    const=not OPT_DEFAULTS.get('test'), default=OPT_DEFAULTS.get('test'),
    help='Run tests.')

  (options, arguments) = parser.parse_args()

  if options.test:
    fail("not available right now")

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
    print "\n"+name+"\n"+str(pos)+": "+cigar
    parse_cigar(cigar, seq, pos)

  if infilehandle is not sys.stdin:
    infilehandle.close()


def parse_cigar(cigar, seq, pos):

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
      print ("insertion after  "+str(ref_pos-1)+": "+seq[alt_pos-4:alt_pos]+"{"
        +seq[alt_pos:alt_pos+length]+"}"+seq[alt_pos+length:alt_pos+length+4])
      alt_pos += length
    elif op == 'D':
      print ("deletion between "+str(ref_pos-1)+" and "+str(ref_pos+length)+": "
        +seq[alt_pos-4:alt_pos]+"["+("-"*length)+"]"
        +seq[alt_pos:alt_pos+4])
      ref_pos += length
    elif op == 'N':
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