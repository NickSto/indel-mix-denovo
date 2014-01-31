#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
import subprocess
import collections
from optparse import OptionParser

OPT_DEFAULTS = {'int':0}
USAGE = """USAGE: %prog [options]
       %prog -mF -f families.tsv -a asm/ -q fastq/ -C -l lastz/ -B bam/asm/"""
DESCRIPTION = """"""
EPILOG = """"""

EXT = {
  'bam':'.bam'
}
# pull out the (full) path base in group 1, extension in group 2
EXT_REGEX = {
  'bam':r'(.+)(\.bam)$',
  'fastq1':r'(.+)(_1\.f(?:ast)?q)$',
  'fastq2':r'(.+)(_2\.f(?:ast)?q)$',
  'fasta':r'(.+)(\.f(?:ast)?a)$',
  'lav':r'(.+)(\.lav)$',
  'tsv':r'(.+)(\.t(sv|ab|xt))$',
}
FILETYPE_EXTS = {
  'fastq1':'fastq1',
  'fastq2':'fastq2',
  'map_fastq1':'fastq1',
  'map_fastq2':'fastq2',
  'bam_in':'bam',
  'bam_out':'bam',
  'asm':'fasta',
  'asm_curated':'fasta',
  'lastz':'lav',
  'report':'tsv',
}
REPORT_KEYS = {
  'h':'contigs-raw',
  'c':'contigs-after',
  'n':'non-ref',
  'd':'duplication',
  'f':'fragmented',
  'F':'failure',
}
REPORT_TYPES = {
  'h':int,
  'c':int,
  'n':int,
  'd':int,
  'f':bool,
  'F':bool,
}

MIN_BWTSW_SIZE = 2000000000 # bytes

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
  parser.add_option('--fastq2-path',
    help="""The second mate reads for paired sequence data. Use when in single-
file mode.""")
  parser.add_option('-Q', '--mapping-fastq-path',
    help="""If you want to use a different set of reads for mapping to the
assembled reference, specify them with this option.""")
  parser.add_option('--mapping-fastq2-path',
    help="""The second mate reads for mapping, when in single-file mode.""")
  parser.add_option('-a', '--asm-path',
    help="""Path to assembly. An assembly is a single FASTA file containing all
the contigs. When the path is a directory, only files ending in .fa or .fasta
will be used.""")
  parser.add_option('-r', '--ref-path',
    help="""Path to reference genome (also a FASTA file). Used to assess and
process the assemblies.""")
  parser.add_option('--curated-asm-path',
    help="""Path to curated assemblies.""")
  parser.add_option('--reports-path',
    help="""Path to assembly reports from asm-curator.py. When path is a
directory, reports must end in .txt, .tsv, or .tab.""")
  parser.add_option('-l', '--lav-path',
    help="""Path to LASTZ output, in LAV format. If provided, the alignment step
will be skipped and this output will be used instead. When the path is a
directory, only files ending in .lav will be used.""")
  parser.add_option('-B', '--bam-path-out',
    help="""Path to file or directory to place alignments to the assembly.""")
  parser.add_option('-C', '--no-curate', action='store_true',
    help="""Do not curate the assemblies before mapping. If family-wise, you
must still provide LAV files to allow determination of which assembly to use.""")
  parser.add_option('-D', '--include-duplications', action='store_true',
    help="""Use assemblies with whole-genome duplications, if necessary.""")
  parser.add_option('-I', '--int', type='int', dest='int',
    default=OPT_DEFAULTS.get('int'),
    help="""default: %default""")
  # parser.set_defaults(verbose=True)

  (options, arguments) = parser.parse_args()

  # if multiple input files, create lists of them
  if options.multi or options.family_wise:
    sample_files = get_multi_files(options)
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

  # Set which pipeline steps are completed, based on input files
  done = get_done_steps(options)

  # If family-wise, loop over families. Else, loop over samples.
  if options.family_wise:
    items = families.keys()
  else:
    items = sample_files.keys()

  # MAIN LOOP
  # Whole pipeline: run once per sample or family
  #TODO: allow samples to fail (expect missing files later)
  for item in items:
    if options.family_wise:
      family = item
      samples = families[family]
    else:
      samples = [item]

    for sample in samples:
      if not done['tofastq']:
        """extract reads from BAM (and dedup)
        input: 'bam_in'
        output: 'fastq1', 'fastq2', 'map_fastq1', 'map_fastq2' """
      if not done['asm']:
        """assemble
        input: 'fastq1', 'fastq2'
        output: 'asm' """
      if not done['lastz']:
        """lastz align to reference
        input: 'asm', ('report'), options.ref_path
        output: 'lastz', ('report') """
        # check options.reports_path or sample_files[sample]['report'] and use
        # that instead of generating another one
      if not done['curated']:
        """curate alignment
        input: 'asm', 'lastz'
        output: 'asm_curated' """

    if not done['align']:
      """choose the best assembly, map to it, merge resulting BAMs
      input: 'asm_curated', 'report', 'map_fastq1', 'map_fastq2'
      output: 'bam_out' """
      # choose best assembly, if multiple
      if options.family_wise:
        best_sample = choose_asm(samples, sample_files, options)
      else:
        best_sample = samples[0]    # ("best_sample" is the *only* sample)
      asmpath = sample_files[best_sample]['asm_curated']
      # align
      for sample in samples:
        fastq1path = sample_files[sample]['map_fastq1']
        fastq2path = sample_files[sample]['map_fastq2']
        bampath = sample_files[sample]['bam_out']
        #TODO: check exit status
        bwa_index(asmpath)
        align_reads(asmpath, fastq1path, fastq2path, bampath)
      # merge, if multiple
      if options.family_wise:
        """merge bam files"""
        #TODO: incorporate family names into sample_files system
        bams = [sample_files[sample]['bam_out'] for sample in samples]
        dirpath = os.path.split(bampath)[0]
        bampath = os.path.join(dirpath, family+'.bam')
        merge_bams(bams, bampath)

    # filter alignment
    # naive variant caller
    # nvc-filter.py
    # inspect-reads.py
    # quick-liftover.py




def get_single_files(options):
  """When not in multi or family mode, pack a dict of single input files.
  See get_multi_files() for structure of dict."""
  #TODO: print errors when directories are given instead of files
  inputs = defaultdict()
  if options.fastq_path and os.path.isfile(options.fastq_path):
    inputs['fastq1'] = options.fastq_path
  if options.fastq2_path and os.path.isfile(options.fastq2_path):
    inputs['fastq2'] = options.fastq2_path
  if options.mapping_fastq_path and os.path.isfile(options.mapping_fastq_path):
    inputs['map_fastq1'] = options.mapping_fastq_path
  if options.mapping_fastq2_path and os.path.isfile(options.mapping_fastq2_path):
    inputs['map_fastq2'] = options.mapping_fastq2_path
  if options.asm_path and os.path.isfile(options.asm_path):
    inputs['asm'] = options.asm_path
  if options.lav_path and os.path.isfile(options.lav_path):
    inputs['lastz'] = options.lav_path
  if options.curated_asm_path and os.path.isfile(options.curated_asm_path):
    inputs['asm_curated'] = options.curated_asm_path
  if options.bam_path_out and os.path.isfile(options.bam_path_out):
    inputs['bam_out'] = options.bam_path_out
  return {'__sample__':inputs}


def get_multi_files(options):
  """Get list of sample names and their associated input files, when inputs are
  given as directories.
  Returns a dict mapping each sample name to a dict of file paths (keyed by the
  type of input)."""
  sample_files = {}
  if options.fastq_path:
    sample_files = get_filetype(sample_files, options.fastq_path, 'fastq1')
    sample_files = get_filetype(sample_files, options.fastq_path, 'fastq2')
    sample_files = get_filetype(sample_files, options.fastq_path, 'map_fastq1')
    sample_files = get_filetype(sample_files, options.fastq_path, 'map_fastq2')
  if options.mapping_fastq_path:
    sample_files = get_filetype(sample_files, options.mapping_fastq_path, 'map_fastq1')
    sample_files = get_filetype(sample_files, options.mapping_fastq_path, 'map_fastq2')
  if options.asm_path:
    sample_files = get_filetype(sample_files, options.asm_path, 'asm')
  if options.lav_path:
    sample_files = get_filetype(sample_files, options.lav_path, 'lastz')
  if options.curated_asm_path:
    sample_files = get_filetype(sample_files, options.curated_asm_path, 'asm_curated')
  if options.reports_path:
    sample_files = get_filetype(sample_files, options.reports_path, 'report')
  if options.bam_path_out:
    sample_files = make_outpaths(sample_files, options.bam_path_out, 'bam_out')
  return sample_files


def get_filetype(sample_files, dirpath, filetype):
  """Get the list of files for a given file type, and add to sample_files dict.
  """
  for filename in os.listdir(dirpath):
    path = os.path.join(dirpath, filename)
    if not os.path.isfile(path):
      continue
    sample = ext_split(filename, FILETYPE_EXTS[filetype])[0]
    if sample:
      files = sample_files.get(sample, {})
      files[filetype] = path
      sample_files[sample] = files
  return sample_files


def make_outpaths(sample_files, dirpath, filetype):
  """"""
  for sample in sample_files:
    extension = EXT[FILETYPE_EXTS[filetype]]
    path = os.path.join(dirpath, sample+extension)
    files = sample_files[sample]
    files[filetype] = path
    sample_files[sample] = files
  return sample_files


def ext_split(path, extype):
  ext_regex = EXT_REGEX[extype]
  match = re.search(ext_regex, path)
  if match:
    return (match.group(1), match.group(2))
  else:
    return (None, None)


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


def get_done_steps(options):
  """"""
  done = collections.defaultdict(bool)
  if (options.fastq_path or options.asm_path or options.lav_path or
      options.curated_asm_path):
    done['tofastq'] = True
  if (options.asm_path or options.lav_path or options.curated_asm_path):
    done['asm'] = True
  if (options.lav_path or options.curated_asm_path):
    done['lastz'] = True
  if (options.curated_asm_path):
    done['curated'] = True
  return done


def choose_asm(samples, sample_files, options):
  """Choose an assembly by prioritizing assembly issues. Some issues are
  dealbreakers, so if all the assemblies have one of these issues, None is
  returned.
  Priority list:
    assembly failure (dealbreaker)
    fragmented assembly (dealbreaker)
    whole-genome duplication (dealbreaker by default)
    non-reference flanks
    number of curated contigs (fewest wins)
    number of original contigs (fewest wins)
  """
  # read in the reports
  reports = {}
  for sample in samples:
    reports[sample] = read_report(sample_files[sample]['report'])
  # narrow the field by prioritized list of assembly issues
  candidates = reports.keys()
  candidates = asms_without('failure', candidates, reports)
  candidates = asms_without('fragmented', candidates, reports)
  if options.include_duplications:
    # if they all have duplication, keep them instead of discarding them
    candidates = narrow_by('duplication', candidates, reports)
  else:
    # eliminate any with duplication, even if it's all of them
    candidates = asms_without('duplication', candidates, reports)
  candidates = narrow_by('non-ref', candidates, reports)
  # of the final candidates, choose the one with the fewest curated contigs
  # (if there's a tie, choose the one with the fewest original contigs)
  candidates.sort(key=lambda sample:
    (reports[sample]['contigs-after'], reports[sample]['contigs-raw'])
  )
  if candidates:
    return candidates[0]
  else:
    return None


def asms_without(issue, samples, reports):
  """Return all samples without the given issue."""
  without = []
  for sample in samples:
    if sample not in reports:
      continue
    if not reports[sample][issue]:
      without.append(sample)
  return without


def narrow_by(issue, samples, reports):
  """Use the given issue to narrow the field, by eliminating assemblies with the
  issue ..UNLESS they all have the issue. Then return the original set of
  assemblies (it's not useful for narrowing the field)."""
  narrowed = asms_without(issue, samples, reports)
  if not narrowed:
    narrowed = samples
  return narrowed


def read_report(reportpath):
  """Read report into a dict, with keys in REPORT_KEYS for each issue, and the
  values given in the report as the values. The dict is a defaultdict, so if an
  issue is not present in the report, it defaults to False."""
  report = collections.defaultdict(bool)
  line_num = 0
  with open(reportpath, 'rU') as reporthandle:
    for line in reporthandle:
      line_num+=1
      line = line.strip()
      if not line:
        continue
      fields = line.split("\t")
      try:
        key = REPORT_KEYS[fields[0]]
        cast_func = REPORT_TYPES[fields[0]]
        report[key] = cast_func(fields[1])
      except (ValueError, KeyError, IndexError):
        fail("Error: invalid report file format on line "+str(line_num)+" of "
          +reportpath)
  return report


def index_ref(ref):
  """Index the reference if it hasn't been already. If indexing was necessary,
  return exit status. Else, return None."""
  indexed = True
  for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']:
    if not os.path.isfile(ref+'.'+ext):
      indexed = False
      break
  if indexed:
    return None
  else:
    print "Indexing assembly "+ref
    command = ['bwa', 'index', '-a']
    if os.path.getsize(ref) > MIN_BWTSW_SIZE:
      command.append('bwtsw')
    else:
      command.append('is')
    command.append(ref)
    return subprocess.call(command)
    print ">>> $ "+' '.join(command)
    

def align_reads(ref, fastq1, fastq2, bam):
  """Align using BWA-MEM, and convert output to indexex and sorted BAM.
  bam parameter is the desired output path.
  If not paired-end data, give None as fastq2."""
  #TODO: check exit statuses
  pathbase = ext_split(bam, 'bam')[0]
  if not pathbase:
    pathbase = bam
  # align
  samtmp = pathbase+'.tmp.sam'
  map_command = ['bwa', 'mem', ref, fastq1]
  if fastq2:
    map_command.append(fastq2)
  with open(samtmp, 'wb') as samhandle:
    status = subprocess.call(map_command, stdout=samhandle)
  print ">>> $ "+' '.join(map_command)+' > '+samtmp
  if status:
    return status
  # convert sam to bam
  bamtmp = pathbase+'.tmp.bam'
  conv_command = ['samtools', 'view', '-Sb', samtmp]
  with open(bamtmp, 'wb') as bamhandle:
    status = subprocess.call(conv_command, stdout=bamhandle)
  print ">>> $ "+' '.join(conv_command)+' > '+bamtmp
  if status:
    return status
  # sort temporary bam
  sort_command = ['samtools', 'sort', bamtmp, pathbase]
  status = subprocess.call(sort_command)
  if status:
    return status
  print ">>> $ "+' '.join(sort_command)
  # index final bam
  index_command = ['samtools', 'index', bam]
  status = subprocess.call(index_command)
  print ">>> $ "+' '.join(index_command)
  # delete temporary bam
  os.remove(samtmp)
  os.remove(bamtmp)
  return status


def merge_bams(bams, merged):
  """Use samtools to merge BAM files.
  merged is the desired output path.
  Reads from different files are labeled with read groups, which are included in
  a valid header."""
  header = make_header(bams)
  command = ['samtools', 'merge', '-r', '-h', header, merged]
  command.extend(bams)
  print '>>> $ '+' '.join(command)
  return subprocess.call(command)


def make_header(bams):
  """Generate a SAM file consisting only of a proper header for a merge of the
  given BAMs.
  The first bam is used as the basis of the header. Read groups in the header
  are taken from filenames in the same way samtools merge does.
  Returns the path to the SAM."""
  # get read group names from BAM filenames, like samtools merge does
  samples = []
  for bam in bams:
    (dirpath, filename) = os.path.split(bam)
    samples.append(ext_split(filename, 'bam')[0])
  header_command = ['samtools', 'view', '-H', bams[0]]
  header = subprocess.check_output(header_command)
  # get list of existing read groups (make sure we don't add a duplicate)
  existing = []
  for line in header.split('\n'):
    if not line.startswith('@RG\t'):
      continue
    fields = line.split('\t')
    for field in fields:
      if field.startswith('ID:'):
        existing.append(field[3:])
  # append new read groups to header
  for sample in samples:
    if sample not in existing:
      header += '@RG\tID:'+sample+'\tSM:'+sample+'\n'
  # print to a temporary file in the same directory as the BAMs
  sampath = os.path.join(dirpath, samples[0]+'.header.sam')
  with open(sampath, 'w') as samhandle:
    samhandle.write(header)
  return sampath


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
