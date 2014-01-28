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

EXT = {
  'bam':'.bam'
}

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
  parser.add_option('-F', '--family-wise', action='store_true',
    help="""Run pipeline on families of samples. Only the best assembly from
each family will be used to map all the samples in that family. The paths to the
input files will be interpreted as directories. If multi mode is not on, all the
files will be assumed to be part of one family.""")
  parser.add_option('-f', '--family-table',
    help="""The path to the family table file. Required when in multi mode.""")
  parser.add_option('-N', '--no-names', dest='has_names', action='store_false',
    default=True,
    help="""Set if the family table does not contain family names as the first
column.""")
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
    sample_files = get_files(options)
  else:
    sample_files = get_single_files(options)

  if options.family_table:
    if not os.path.isfile(options.family_table):
      fail("Error: cannot open family table "+options.family_table)
    families = read_family_table(options.family_table, options.has_names)
  elif options.family_wise and options.multi:
    parser.print_help()
    fail("\nError: If running on multiple families, you must provide a family "
      +"table.")
  else:
    families = {'__family__':['__sample__']}

  # If family-wise, loop over families. Else, loop over samples.
  if options.family_wise:
    items = families.keys()
  else:
    items = sample_files.keys()

  for item in items:
    if options.family_wise:
      family = item
      samples = families[family]
    else:
      samples = [item]

    for sample in samples:
      fastq1path = sample_files[sample]['fastq1']
      fastq2path = sample_files[sample]['fastq2']
      asmpath = sample_files[sample]['asm']
      # extract reads from BAM
      # assemble
      # lastz align to reference
      # curate alignment

    if options.family_wise:
      best_sample = choose_asm(samples, sample_files)
    else:
      best_sample = samples[0]

    asmpath = sample_files[best_sample]['asm']
    for sample in samples:
      fastq1path = sample_files[sample]['fastq1']
      fastq2path = sample_files[sample]['fastq2']
      map_to_asm(asmpath, fastq1path, fastq2path, bam)

    if options.family_wise:
      """merge bam files"""
    else:
      pass

    # filter alignment
    # naive variant caller
    # nvc-filter.py
    # inspect-reads.py
    # quick-liftover.py


    print sample+':'
    for filetype in ['asm', 'fastq1', 'fastq2', 'lav']:
      print filetype+":\t"+inputs.get(filetype)
    # map fastq reads to asm fasta



def get_single_files(options):
  """When not in multi or family mode, pack a dict of single input files.
  See get_files() for structure of dict."""
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
  return {'__sample__':inputs}


def get_files(options):
  """Get list of sample names and their associated input files.
  Returns a dict mapping each sample name to a dict of file paths (keyed by the
  type of input)."""
  sample_files = {}
  if options.asm_path:
    sample_files = get_filetype(sample_files, options.asm_path, 'asm')
  if options.fastq_path:
    sample_files = get_filetype(sample_files, options.fastq_path, 'fastq1')
    sample_files = get_filetype(sample_files, options.fastq_path, 'fastq2')
  if options.asm_path:
    sample_files = get_filetype(sample_files, options.lav_path, 'lav')
  if options.bam_path_out:
    sample_files = make_outpaths(sample_files, options.bam_path_out, 'bam')
  return sample_files


def get_filetype(sample_files, dirpath, filetype):
  """Get the list of files for a given file type, and add to sample_files dict.
  """
  for filename in os.listdir(dirpath):
    if filetype.startswith('fastq'):
      print "file in "+dirpath+": "+filename
    path = os.path.join(dirpath, filename)
    if not os.path.isfile(path):
      continue
    match = re.search(EXT_REGEX[filetype], filename)
    if match:
      sample = match.group(1)
      files = sample_files.get(sample, {})
      files[filetype] = path
      sample_files[sample] = files
  return sample_files


def make_outpaths(sample_files, dirpath, filetype):
  """"""
  for sample in sample_files:
    path = os.path.join(dirpath, sample+EXT[filetype])
    files = sample_files[sample]
    files[filetype] = path
    sample_files[sample] = files
  return sample_files


def read_family_table(tablepath, names=True):
  families = {}
  count = 0
  with open(tablepath, 'rU') as tablehandle:
    for line in tablehandle:
      samples = []
      line = line.strip()
      if not line:
        continue
      columns = line.split('\t')
      if names:
        name = columns[0]
        columns = columns[1:]
      else:
        name = str(count)
        count+=1
      for sample in columns:
        samples.append(sample)
      families[name] = samples
  return families


def map_to_asm(asm, fastq1, fastq2, bam):
  """"""

def choose_asm(samples, sample_files):
  """"""
  return samples[0]


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
