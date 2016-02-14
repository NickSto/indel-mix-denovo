#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import argparse
import pipeline

ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('-1', '--fastqs1', metavar='path/to/sampleA_1.fq,path/to/sampleB_1.fq',
    help='The input FASTQ files, mate 1. Give a comma-separated list of paths to the files.')
  parser.add_argument('-2', '--fastqs2', metavar='path/to/sampleA_2.fq,path/to/sampleB_2.fq',
    help='The input FASTQ files, mate 2. Give a comma-separated list of paths to the files.')
  parser.add_argument('-s', '--samples', metavar='sampleA,sampleB,sampleC',
    help='The sample names, comma-separated.')
  parser.add_argument('-f', '--family',
    help='The family name.')
  parser.add_argument('-i', '--input-dir', metavar='path/to/input/directory')
  parser.add_argument('-L', '--family-file', metavar='path/to/families.tsv')

  parser, opts, param = pipeline.add_parameters(parser)
  opts_dict = opts_to_dict(opts)
  args = parser.parse_args(argv[1:])

  for arg in filter(lambda arg: not arg.startswith('_'), dir(args)):
    if arg not in opts_dict:
      continue
    opt = opts_dict[arg]
    arg_value = getattr(args, arg)
    if arg_value != opt.default:
      #TODO: Handle boolean flags.
      opt_strings.append(opt.option_strings[0])
      opt_strings.append(str(arg_value))


def opts_to_dict(opts):
  opts_dict = {}
  for opt in opts:
    opts_dict[opt.dest] = opt
  return opts_dict


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
