#!/usr/bin/env python
import re
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'test':False}
USAGE = "USAGE: %prog [options]"
DESCRIPTION = """"""
EPILOG = """"""

TESTS = [
  '251M',       # M
  '3S248M',     # SM
  '1M2I248M',   # MIM
  '166M1D85M',  # MDM
  '179M62I2M2D2M6S',  # MIMDMS
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
    cigars = TESTS
  else:
    cigars = arguments

  if not cigars:
    parser.print_help()
    fail("Provide one or more CIGAR strings as arguments")

  for cigar in cigars:
    print "\n"+cigar
    ref_pos = 0
    alt_pos = 0
    for move in re.findall(r'[^A-Z]*[A-Z]', cigar):
      print move
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
          fail("Error: CIGAR contains movement in unrecognized format.\n"
            +"CIGAR: "+cigar+"\tmovement: "+move)
      if op == 'M':
        ref_pos += length
        alt_pos += length
      elif op == 'I':
        print "insertion at "+str(ref_pos)
        alt_pos += length
      elif op == 'D':
        print "deletion from "+str(ref_pos)+" to "+str(ref_pos+length)
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