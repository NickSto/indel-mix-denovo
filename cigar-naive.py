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
    print cigar
    pos = 0
    for move in re.findall(r'[^A-Z]*[A-Z]', cigar):
      match = re.search(r'^(\d+)([MIDSHPN])$', move)
      if match:
        length = match.group(1)
        type = match.group(2)
      else:
        match = re.search(r'^([MIDSHPN])$', move)
        if match:
          length = 1
          type = match.group(1)
        else:
          fail("Error: CIGAR contains movement in unrecognized format.\n"
            +"CIGAR: "+cigar+"\tmovement: "+move)
      print type+" "+str(length)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()