#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import argparse

ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  opts = []
  opts.append(parser.add_argument('-1', '--fastqs1',
    metavar='path/to/sampleA_1.fq,path/to/sampleB_1.fq',
    help='The input FASTQ files, mate 1. Give a comma-separated list of paths to the files.'))
  opts.append(parser.add_argument('-2', '--fastqs2',
    metavar='path/to/sampleA_2.fq,path/to/sampleB_2.fq',
    help='The input FASTQ files, mate 2. Give a comma-separated list of paths to the files.'))
  opts.append(parser.add_argument('-s', '--samples', metavar='sampleA,sampleB,sampleC',
    help='The sample names, comma-separated.'))
  opts.append(parser.add_argument('-f', '--family',
    help='The family name.'))
  opts.append(parser.add_argument('-i', '--input-dir', metavar='path/to/input/directory'))
  opts.append(parser.add_argument('-L', '--family-file', metavar='path/to/families.tsv'))

  opts_dict = opts_to_dict(opts)
  our_argv, pipeline_argv = separate_argv(argv, opts_dict)
  args = parser.parse_args(our_argv)


def opts_to_dict(opts):
  opts_dict = {}
  for opt in opts:
    for option_string in opt.option_strings:
      opts_dict[option_string] = opt
  return opts_dict


def separate_argv(argv, opts_dict):
  our_argv = []
  pipeline_argv = []
  nargs = 0
  for arg in argv[1:]:
    if arg in opts_dict:
      our_argv.append(arg)
      opt = opts_dict[arg]
      if opt.nargs is None:
        nargs = 1
      else:
        nargs = opt.nargs
    elif nargs > 0:
      our_argv.append(arg)
      nargs -= 1
    elif arg in ('-h', '--help'):
      our_argv.append(arg)
    else:
      pipeline_argv.append(arg)
  return our_argv, pipeline_argv


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
