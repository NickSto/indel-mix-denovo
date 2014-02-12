#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import itertools

OPT_DEFAULTS = {'ending':'filt'}
USAGE = "%(prog)s [options] -e ending infile1.tsv [infile2.tsv [..]]"
DESCRIPTION = """Filter variant files as a group. E.g. only remove a site if it
fails a filter in *all* files. *IMPORTANT* All files must contain the same sites
on the same lines. If site locations are not the same for a given line number,
this will fail."""
EPILOG = """"""

def main():

  parser = argparse.ArgumentParser(
    description=DESCRIPTION, usage=USAGE, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infiles', metavar='infile.tsv', nargs='+',
    help="""Input file. Format must be the tsv output of inspect-reads.py.""")
  parser.add_argument('-e', '--ending',
    help="""The ending to append to the output filename (before the extension.
      E.g. "-e nobias" on file "R19S1.tsv" will write to the file
      "R19S1-nobias.tsv". Default: "%(default)s\"""")
  parser.add_argument('-s', '--strand-bias', metavar='BIAS', type=float,
    help="""Strand bias threshold. Sites with a bias above this value fail.
      The bias statistic is method 1 (SB) of Guo et al., 2012.""")
  parser.add_argument('-m', '--mate-bias', metavar='BIAS', type=float,
    help="""Mate bias threshold. Sites with a bias above this value fail.
      The bias statistic is calculated identically to the strand bias, replacing
      forward/reverse strands with first/second mate in the pair.""")
  parser.add_argument('-a', '--any', action='store_true',
    help="""***NOT YET IMPLEMENTED*** Remove a site if fails a filter in *any*
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
  done = False
  linenum = 0
  while not done:
    linenum+=1
    last_site = None
    passed = False
    # check each infile to see if the site passes the filters
    lines = []
    for infile in infiles:
      place = (linenum, infile.name)
      line = infile.readline()
      lines.append(line)
      # EOF
      if not line:
        done = True
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
        passed = passes(fields, args)
      except ValueError:
        format_fail("Wrong type found when parsing line %d of file %s.", place)
      except IndexError:
        format_fail("Too few columns on line %d of file %s.", place)
      # sys.stdout.write(infile.name+'\t')
    if done:
      break
    if passed:
      for (line, outfile) in zip(lines, outfiles):
        outfile.write(line)
    # print "%2d: %s" % (linenum, fields[1])

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
