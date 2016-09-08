#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import argparse

ARG_DEFAULTS = {'input':sys.stdin, 'method':'range', 'thres':3}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='vars_asm.tsv', type=argparse.FileType('r'), nargs='?',
    help='Output of inspect-reads.py')
  parser.add_argument('-m', '--method', choices=('range', 'r', 'rsquared'))
  parser.add_argument('-f', '--factor', type=float,
    help='Default: %(default)s')
  parser.add_argument('-r', '--max-r', type=float,
    help='Default: %(default)s')

  args = parser.parse_args(argv[1:])

  for line in args.input:
    if line.startswith('#'):
      continue
    fields = line.rstrip('\r\n').split('\t')
    if fields[0] == 'Sample':
      continue
    pos_dist_str = fields[22]
    total = int(fields[6])
    pos_dist = [int(n) for n in pos_dist_str.split(',')]
    bins = len(pos_dist)
    expected = total/bins
    min_count = expected/args.thres
    max_count = args.thres*expected
    r = 0
    failed = False
    for count in pos_dist:
      if args.method == 'range':
        if count < min_count or count > max_count:
          failed = True
          break
      elif args.method == 'r':
        r += abs(count - expected)
    if args.method == 'r':
      if r/total > args.thres:
        failed = True
    if failed:
      pos_dist_ratio = ' '.join(['{:4.1f}'.format(n/expected) for n in pos_dist])
      pos_dist_raw   = ' '.join(['{:3d}'.format(n) for n in pos_dist])
      print('{2:9s}{3:5s}{4:10s}{5:9s}{7:9s}{19:10s}{20:10s}{6:5s}|  {raw:43s}'
            .format(*fields, ratio=pos_dist_ratio, raw=pos_dist_raw))


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
