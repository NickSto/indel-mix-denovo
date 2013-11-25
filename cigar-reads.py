#!/usr/bin/env python
import re
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'bool':False}
USAGE = "USAGE: %prog reads.sam cigarpattern [cigarpattern [cigarpattern [..]]]"
DESCRIPTION = """Get alignments (reads) which have a given cigar pattern"""
EPILOG = """"""

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-s', '--str', dest='str',
    default=OPT_DEFAULTS.get('str'), help='default: %default')
  parser.add_option('-i', '--int', dest='int', type='int',
    default=OPT_DEFAULTS.get('int'), help='')
  parser.add_option('-f', '--float', dest='float', type='float',
    default=OPT_DEFAULTS.get('float'), help='')
  parser.add_option('-b', '--bool', dest='bool', action='store_const',
    const=not OPT_DEFAULTS.get('bool'), default=OPT_DEFAULTS.get('bool'),
    help='')

  (options, arguments) = parser.parse_args()

  if not arguments:
    parser.print_help()
    fail("Give a CIGAR pattern")
  else:
    samfile = arguments[0]
    cigars = arguments[1:]

  pattern_regexes = []
  for cigar in cigars:
    pattern = re.sub(r'[^MIDNSHPX=]', '', cigar)
    if not pattern:
      fail("Error: Invalid CIGAR pattern: "+str(cigar))
    pattern_regex = r'^'
    for op in pattern:
      pattern_regex += r'\d*'+op
    pattern_regex += r'$'
    pattern_recomp = re.compile(pattern_regex)
    pattern_regexes.append(pattern_recomp)

  if samfile == '-':
    samhandle = sys.stdin
  else:
    samhandle = open(samfile, 'rU')

  for line in samhandle:
    if line[0] == '@' and re.search(r'^@[A-Z]{2}\t[A-Z]{2}', line):
      continue
    fields = line.split("\t")
    for pattern in pattern_regexes:
      if pattern.match(fields[5]):
        sys.stdout.write(line)
        break

  if isinstance(samhandle, file):
    samhandle.close()


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()