#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import itertools

OPT_DEFAULTS = {'ending':'filt'}
USAGE = "%(prog)s [options] -e ending infile1.tsv [infile2.tsv [..]]"
DESCRIPTION = """Filter variant files as a group. By default, if a site passes
all filters in at least one file, it will be kept (in all files)."""
EPILOG = """*IMPORTANT* All files must contain the same sites
on the same lines. If site locations are not the same for a given line number,
this will fail."""

def main():

  parser = argparse.ArgumentParser(
    description=DESCRIPTION, usage=USAGE, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infiles', metavar='infile.tsv', nargs='+',
    help="""Input file. Format must be the tsv output of inspect-reads.py.""")
  parser.add_argument('-e', '--ending',
    help="""The ending to append to the output filename (before the extension.
      E.g. "-e nobias" on file "R19S1.tsv" will write to the file
      "R19S1-nobias.tsv". Default: "%(default)s".""")
  parser.add_argument('-s', '--strand-bias', metavar='BIAS', type=float,
    help="""Strand bias threshold. Sites with a bias above this value fail.
      The bias statistic is method 1 (SB) of Guo et al., 2012.""")
  parser.add_argument('-m', '--mate-bias', metavar='BIAS', type=float,
    help="""Mate bias threshold. Sites with a bias above this value fail.
      The bias statistic is calculated identically to the strand bias, replacing
      forward/reverse strands with first/second mate in the pair.""")
  parser.add_argument('-a', '--any', action='store_true',
    help="""**Not tested yet** Remove a site if fails a filter in *any*
      file (default is to remove if it fails in *all* files).""")

  args = parser.parse_args()

  # set as global to allow access from cleanup() function
  global infiles
  global outfiles
  global success
  infiles = []
  outfiles = []
  success = False

  # open input and output filehandles
  for infile in args.infiles:
    infiles.append(open(infile, 'rU'))
    outfile = get_outfile_name(infile, args)
    outfiles.append(open(outfile, 'w'))

  # main loop
  eof = False
  linenum = 0
  while not eof:
    linenum+=1
    last_site = None
    if args.any:
      all_passed = True
    else:
      all_passed = False
    # check each infile to see if the site passes the filters
    lines = []
    for infile in infiles:
      place = (linenum, infile.name)
      line = infile.readline()
      lines.append(line)
      if not line:
        eof = True
        break
      fields = line.strip().split('\t')
      # format check; only require columns that are actually needed
      if len(fields) < 2:
        format_fail("Too few columns on line %d of file %s.", place)
      # site matches up with the ones in the other files?
      this_site = (fields[0], fields[1])
      if last_site and last_site != this_site:
        fail("Error: lines do not match up (mismatching site on line %d of "
          "file %s)." % place)
      last_site = this_site
      # apply filters
      try:
        this_passed = passes(fields, args)
      except ValueError:
        format_fail("Wrong type found when parsing line %d of file %s.", place)
      except IndexError:
        format_fail("Too few columns on line %d of file %s.", place)
      if args.any:
        all_passed = all_passed and this_passed
      else:
        all_passed = all_passed or this_passed
    if eof:
      break
    if all_passed:
      for (line, outfile) in zip(lines, outfiles):
        outfile.write(line)

  success = True
  cleanup()


def get_outfile_name(infile, args):
  (root, ext) = os.path.splitext(infile)
  base = root+'-'+args.ending
  count = 0
  name = base+ext
  while os.path.exists(name):
    count+=1
    name = base+str(count)+ext
    if count > 1000:
      fail("Error: How many files can you have named "+base+ext+"?")
  return name


def passes(fields, args):
  passed = True
  if args.strand_bias:
    if float(fields[18]) > args.strand_bias:
      passed = False
  if args.mate_bias:
    if float(fields[19]) > args.mate_bias:
      passed = False
  return passed


def cleanup():
  for infile in infiles:
    if isinstance(infile, file):
      infile.close()
  for outfile in outfiles:
    if isinstance(outfile, file):
      outfile.close()
  # delete output files if it aborted mid-way
  if not success:
    for outfile in outfiles:
      os.remove(outfile.name)


def format_fail(message, place):
  fail_message = "Format Error: "+message % place
  fail(fail_message)


def fail(message):
  cleanup()
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  main()
