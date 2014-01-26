#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'int':0}
USAGE = """USAGE: %prog [options]
       %prog -mF -f families.tsv -a asm/ -q fastq/ -C -l lastz/ -B bam/asm/"""
DESCRIPTION = """"""
EPILOG = """"""

# pull out the filename base in group 1, extension in group 2
EXT_REGEX = {
  'fastq':r'([^/]+)(_[12]\.f(?:ast)?q)$',
  'fasta':r'([^/]+)(\.f(?:ast)?a)$',
}

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-m', '--multi', action='store_true',
    help="""Run the pipeline on multiple samples. When in multi mode, the paths
to the input files (except the family table) must be directories containing the
files. Files from the same sample must have identical base filenames (everything
except the extension) in each directory.""")
  parser.add_option('-F', '--family-wise',
    help="""Run pipeline on families of samples. Only the best assembly from
each family will be used to map all the samples in that family. The paths to the
input files will be interpreted as directories. If multi mode is not on, all the
files will be assumed to be part of one family.""")
  parser.add_option('-f', '--family-table',
    help="""The path to the family table file. Required when in multi mode.""")
  parser.add_option('-o', '--output-path',
    help="""Path to output file or directory.""")
  parser.add_option('-q', '--fastq-path',
    help="""Path to input sequence data. When the path is a directory, only
files ending in .fq or .fastq will be used. Paired-end filenames must end in
_1.fq and _2.fq (or .fastq) to designate pairs.""")
  parser.add_option('-a', '--asm-path',
    help="""Path to assembly. An assembly is a single FASTA file containing all
the contigs. When the path is a directory, only files ending in .fa or .fasta
will be used.""")
  parser.add_option('-B', '--bam-path-out',
    help="""Path to file or directory to place alignments to the assembly.""")
  parser.add_option('-C', '--no-curate', action='store_true',
    help="""Do not curate the assemblies before mapping. If family-wise, you
must still provide LAV files to allow determination of which assembly to use.""")
  parser.add_option('-l', '--lav-path',
    help="""Path to LASTZ output, in LAV format. If provided, the alignment step
will be skipped and this output will be used instead. When the path is a
directory, only files ending in .lav will be used.""")
  parser.add_option('-I', '--int', type='int', dest='int',
    default=OPT_DEFAULTS.get('int'),
    help="""default: %default""")
  # parser.set_defaults(verbose=True)

  (options, arguments) = parser.parse_args()

  # if multiple input files, create lists of them
  if options.multi or options.family_wise:
    if options.fastq_path:
      basenames = get_basenames(options.fastq_path, 'fastq')
    elif options.asm_path:
      basenames = get_basenames(options.asm_path, 'fasta')
    else:
      parser.print_help()
      fail("\nError: You must provide correct input datasets.")
  else:
    fail("Error: Single-file mode not yet implemented.")

  #TODO: check that all the filenames correspond



def get_basenames(dirpath, filetype):
  """Return basenames of files in dirpath in a set."""
  basenames = set()
  for filename in os.listdir(dirpath):
    if not os.path.isfile(os.path.join(dirpath, filename)):
      continue
    match = re.search(EXT_REGEX[filetype], filename)
    if match:
      basenames.add(match.group(1))
  return basenames


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
