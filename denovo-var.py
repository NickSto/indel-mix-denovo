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
  'bam':r'([^/]+)(\.bam)$',
  'fastq1':r'([^/]+)(_1\.f(?:ast)?q)$',
  'fastq2':r'([^/]+)(_2\.f(?:ast)?q)$',
  'asm':r'([^/]+)(\.f(?:ast)?a)$',
  'lav':r'([^/]+)(\.lav)$',
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
  parser.add_option('-Q', '--fastq2-path',
    help="""The second mate reads for paired sequence data. Use when in single-
file mode.""")
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
    input_files = get_input_files(options)
  else:
    input_files = get_single_input_files(options)

  for (basename, inputs) in input_files.items():
    print basename+':'
    for filetype in ['asm', 'fastq1', 'fastq2', 'lav']:
      print filetype+":\t"+inputs.get(filetype)

  #TODO: check that all the filenames correspond



def get_single_input_files(options):
  """When not in multi or family mode, pack a dict of single input files.
  See get_input_files() for structure of dict."""
  #TODO: print errors when directories are given instead of files
  inputs = {}
  if options.asm_path and os.path.isfile(options.asm_path):
    inputs['asm'] = options.asm_path
  if options.fastq_path and os.path.isfile(options.fastq_path):
    inputs['fastq1'] = options.fastq_path
  if options.fastq2_path and os.path.isfile(options.fastq2_path):
    inputs['fastq2'] = options.fastq2_path
  if options.lav_path and os.path.isfile(options.lav_path):
    inputs['lav'] = options.lav_path
  return {'__base__':inputs}


def get_input_files(options):
  """Get list of file basenames and their associated input files.
  Returns a dict mapping each basename to a dict of file paths (keyed by the
  type of input)."""
  input_files = {}
  if options.asm_path:
    input_files = get_type_of_files(input_files, options.asm_path, 'asm')
  if options.fastq_path:
    print "gathering fastq's from "+options.fastq_path
    input_files = get_type_of_files(input_files, options.fastq_path, 'fastq1')
    input_files = get_type_of_files(input_files, options.fastq_path, 'fastq2')
  if options.asm_path:
    input_files = get_type_of_files(input_files, options.lav_path, 'lav')
  return input_files


def get_type_of_files(input_files, dirpath, filetype):
  """Get the list of files for a given type of input, and add to input_files
  dict."""
  for filename in os.listdir(dirpath):
    if filetype.startswith('fastq'):
      print "file in "+dirpath+": "+filename
    path = os.path.join(dirpath, filename)
    if not os.path.isfile(path):
      continue
    match = re.search(EXT_REGEX[filetype], filename)
    if match:
      basename = match.group(1)
      inputs = input_files.get(basename, {})
      inputs[filetype] = path
      input_files[basename] = inputs
      if filetype.startswith('fastq'):
        print "added "+filename+" to "+basename
  return input_files


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
