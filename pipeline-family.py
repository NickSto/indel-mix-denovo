#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import logging
import argparse

FASTQ_EXTS = ('.fq', '.fastq')

ARG_DEFAULTS = {'log':sys.stderr}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  opts = []
  parser.add_argument('-r', '--ref', metavar='reference.fa',
    help='The reference genome.')
  parser.add_argument('-o', 'outdir', metavar='output/directory/path',
    help='Destination directory to place the output. If running on a single family, this will '
         'create a separate subdirectory for each sample in this output directory. If running on '
         'multiple families, it will create a separate subdirectory for each family, each '
         'containing a subdirectory for each sample in the family.')
  opts.append(parser.add_argument('-1', '--fastqs1',
    metavar='path/to/sampleA_1.fq,path/to/sampleB_1.fq',
    help='The input FASTQ files, mate 1. Give a comma-separated list of paths to the files. '
         'Optional. If not provided, the paths will be inferred using --input-dir and --samples, '
         'assuming the files are named [sample]_1.fq, etc.'))
  opts.append(parser.add_argument('-2', '--fastqs2',
    metavar='path/to/sampleA_2.fq,path/to/sampleB_2.fq',
    help='The input FASTQ files, mate 2. Give a comma-separated list of paths to the files.'))
  opts.append(parser.add_argument('-s', '--samples', metavar='sampleA,sampleB,sampleC',
    help='The sample names, comma-separated.'))
  opts.append(parser.add_argument('-f', '--family',
    help='The family name.'))
  opts.append(parser.add_argument('-i', '--input-dir', metavar='path/to/input/directory',
    help='The directory where the FASTQ files are stored. Only needed if --fastqs1 and --fastqs2 '
         'are not provided.'))
  opts.append(parser.add_argument('-L', '--family-file', metavar='path/to/families.tsv'))
  opts.append(parser.add_argument('-q', '--quiet', dest='log_level', action='store_const',
    const=logging.ERROR,
    help='Print messages only on terminal errors.'))
  opts.append(parser.add_argument('-v', '--verbose', dest='log_level', action='store_const',
    const=logging.INFO,
    help='Print informational messages in addition to warnings and errors.'))
  opts.append(parser.add_argument('-D', '--debug', dest='log_level', action='store_const',
    const=logging.DEBUG,
    help='Turn debug messages on.'))
  opts.append(parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Will overwrite if it exists.'))

  opts_dict = opts_to_dict(opts)
  our_argv, pipeline_argv = separate_argv(argv, opts_dict)
  args = parser.parse_args(our_argv)

  logging.basicConfig(stream=args.log, level=args.log_level, format='%(levelname)s: %(message)s')
  tone_down_logger()

  script_dir = os.path.dirname(os.path.realpath(__file__))

  fastqs = get_fastqs(args.fastqs1, args.fastqs2, args.samples, args.input_dir, FASTQ_EXTS)
  if fastqs is None:
    parser.print_help()
    return 1

  # {script_dir}/pipeline.py {args.ref} {fastq1} {fastq2} {outdir} -s {sample} {' '.join(pipeline_argv)}


def opts_to_dict(opts):
  """Convert a list of argparse._StoreAction's to a dict mapping each option_string to the _StoreAction
  it's associated with."""
  opts_dict = {}
  for opt in opts:
    for option_string in opt.option_strings:
      opts_dict[option_string] = opt
  return opts_dict


def separate_argv(argv, opts_dict):
  """Parse the raw arguments list and separate the ones for this script from the ones for pipeline.py.
  Returns the arguments in two separate lists, respectively."""
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


def get_fastqs(args_fastqs1, args_fastqs2, args_samples, args_input_dir, fastq_exts):
  """Parse input arguments and return a list of fastq paths.
  Gets the paths either directly from --fastqs1 and --fastqs2 or by combining the --input-dir with
  the --samples, plus "_1.fq" and "_2.fq". Also checks that the paths exist.
  Returns a list of 2-tuples. Each tuple is a matched pair of fastq paths."""
  if args_fastqs1 and args_fastqs2:
    fastqs1 = args_fastqs1.split(',')
    fastqs2 = args_fastqs2.split(',')
    for fastq in fastqs1 + fastqs2:
      if not os.path.isfile(fastq):
        raise IOError('Could not find input FASTQ file "'+fastq+'".')
  elif args_samples and args_input_dir:
    # Fall back to --samples and --input-dir if either --fastqs1 or --fastqs2 is omitted.
    samples = args_samples.split(',')
    fastqs1 = []
    fastqs2 = []
    for sample in samples:
      fastq1_base = os.path.join(args_input_dir, sample+'_1')
      fastqs1.append(find_fastq(fastq1_base, fastq_exts))
      fastq2_base = os.path.join(args_input_dir, sample+'_2')
      fastqs2.append(find_fastq(fastq2_base, fastq_exts))
  else:
    logging.error('Must provide either --fastqs1 and --fastqs2 or --input-dir and --samples.')
    return None
  return zip(fastqs1, fastqs2)


def find_fastq(fastq_base, fastq_exts):
  """Return the path to a fastq file formed by fastq_base + fastq_ext.
  Will loop through all fastq_exts until it finds a file that exists. If none succeeds, it will
  raise an IOError."""
  for fastq_ext in fastq_exts:
    fastq = fastq_base + fastq_ext
    if os.path.isfile(fastq):
      return fastq
  raise IOError('Could not find input FASTQ file "'+fastq_base+fastq_exts[0]+'". When using '
                '--input-dir and --samples, make sure the FASTQ paths follow the format [input_dir]'
                '/[sample]_1'+fastq_exts[0]+'.')


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
