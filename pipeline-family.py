#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import logging
import argparse
import subprocess
import collections

FASTQ_EXTS = ('.fq', '.fastq')

ARG_DEFAULTS = {'log':sys.stderr, 'max_contigs':500}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""

REPORT_CODES = {
  'h':'contigs-raw',
  'c':'contigs-after',
  'g':'gaps',
  'n':'non-ref',
  'd':'duplication',
  'f':'fragmented',
  'F':'failure',
}


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  opts = []
  opts.append(parser.add_argument('-r', '--ref', metavar='reference.fa', required=True,
    help='The reference genome.'))
  opts.append(parser.add_argument('-o', '--outdir', metavar='output/directory/path', required=True,
    help='Destination directory to place the output. This will create a separate subdirectory for '
         'each sample in this output directory. The subdirectories must either be empty or not '
         'exist yet.'))
         # 'If running on multiple families, it will create a separate subdirectory for each '
         # 'family, each containing a subdirectory for each sample in the family.')
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
  opts.append(parser.add_argument('-b', '--to-script', metavar='script-name.sh',
    help='Write the commands to scripts with this name in each sample directory.'))
  opts.append(parser.add_argument('-C', '--max-contigs',
    help='The maximum allowed number of contigs per assembly (approximated by the number of LASTZ '
         'hits to the reference). Set to 0 to allow any number. Default: %(default)s'))
  opts.append(parser.add_argument('-D', '--include-dup', action='store_true',
    help='Don\'t fail if all assemblies in the family show whole-genome duplications.'))
  # opts.append(parser.add_argument('-L', '--family-file', metavar='path/to/families.tsv'))
  # Logging settings
  opts.append(parser.add_argument('-q', '--quiet', dest='log_level', action='store_const',
    const=logging.ERROR,
    help='Print messages only on terminal errors.'))
  opts.append(parser.add_argument('-v', '--verbose', dest='log_level', action='store_const',
    const=logging.INFO,
    help='Print informational messages in addition to warnings and errors.'))
  opts.append(parser.add_argument('--debug', dest='log_level', action='store_const',
    const=logging.DEBUG,
    help='Turn debug messages on.'))
  opts.append(parser.add_argument('--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Will overwrite if it exists.'))

  opts_dict = opts_to_dict(opts)
  our_args, pipeline_args = separate_argv(argv, opts_dict)
  args = parser.parse_args(our_args)

  logging.basicConfig(stream=args.log, level=args.log_level, format='%(levelname)s: %(message)s')
  tone_down_logger()

  script_dir = os.path.dirname(os.path.realpath(__file__))

  fastqs1, fastqs2, samples = get_fastqs(args.fastqs1, args.fastqs2, args.samples, args.input_dir,
                                         FASTQ_EXTS)
  if fastqs1 is None:
    parser.print_help()
    return 1

  # Check and create output directories.
  outdirs = []
  for sample in samples:
    outdir = os.path.join(args.outdir, sample)
    outdirs.append(outdir)
    if os.path.exists(outdir) and not os.path.isdir(outdir):
      logging.error('Output directory path "{}" already exists but is not a directory.'
                    .format(outdir))
      return 1
    if os.path.isdir(outdir):
      if os.listdir(outdir):
        logging.error('Output directory path "{}" already exists but is not empty.'.format(outdir))
        return 1
    else:
      os.makedirs(outdir)

  # Check and overwrite the --to-script.
  if args.to_script:
    if os.path.split(args.to_script)[0] != '':
      logging.error('--to-script should be a filename, not a path ({}).')
      return 1

  # Run first half of pipeline.
  run_pipelines(script_dir, args.ref, samples, fastqs1, fastqs2, outdirs, pipeline_args, end=3,
                script=args.to_script)

  # Check the assemblies and choose the best.
  reports = {}
  for sample, outdir in zip(samples, outdirs):
    lav_raw = os.path.join(outdir, 'lav_raw.lav')
    command = 'lastz {ref} {outdir}/asmdir/contigs.fasta > {lav_raw}'.format(ref=args.ref,
                                                                             outdir=outdir,
                                                                             lav_raw=lav_raw)
    subprocess.call(command, shell=True)
    command = '{script_dir}/asm-curator.py -l {lav_raw}'.format(script_dir=script_dir,
                                                                lav_raw=lav_raw)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    reports[sample] = parse_report(process.communicate()[0], REPORT_CODES)
  best_sample = choose_asm(reports, args.max_contigs, args.include_dup)
  if best_sample is None:
    logging.error('No high-quality assembly chosen.')
    return 1
  best_asm = os.path.join(args.outdir, best_sample, 'asm.fa')
  best_lav = os.path.join(args.outdir, best_sample, 'lav.lav')
  # Replace the assembly/lav files in each sample directory with a link to the best one.
  for outdir in outdirs:
    asm = os.path.join(outdir, 'asm.fa')
    lav = os.path.join(outdir, 'lav.lav')
    if asm == best_asm and lav == best_lav:
      continue
    if os.path.exists(asm):
      os.rename(asm, os.path.join(outdir, 'asm.orig.fa'))
    if os.path.exists(lav):
      os.rename(lav, os.path.join(outdir, 'lav.orig.lav'))
    relative_link(asm, best_asm)
    relative_link(lav, best_lav)

  # Align the reads to the assembly.
  # This is done separately from the rest of the 2nd half of the pipeline because we need to give
  # different values as the "sample name" in the two cases. Here, it's the actual sample names.
  # Later, it's the name of the chosen assembly.
  run_pipelines(script_dir, args.ref, samples, fastqs1, fastqs2, outdirs, pipeline_args, begin=4,
                end=4, script=args.to_script)

  # Run second half of pipeline.
  best_samples = (best_sample,) * len(samples)
  run_pipelines(script_dir, args.ref, best_samples, fastqs1, fastqs2, outdirs, pipeline_args,
                begin=5, script=args.to_script)


def opts_to_dict(opts):
  """Convert a list of argparse._StoreAction's to a dict mapping each option_string
  to the _StoreAction it's associated with."""
  opts_dict = {}
  for opt in opts:
    for option_string in opt.option_strings:
      opts_dict[option_string] = opt
  return opts_dict


def separate_argv(argv, opts_dict):
  """Parse the raw arguments list and separate the ones for this script from the ones
  for pipeline.py. Returns the arguments in two separate lists, respectively."""
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
  Returns a list of 3-tuples. Each tuple contains (fastq1 path, fastq2 path, sample name)."""
  if args_fastqs1 and args_fastqs2:
    fastqs1 = args_fastqs1.split(',')
    fastqs2 = args_fastqs2.split(',')
    # Check that the files exist.
    for fastq in fastqs1 + fastqs2:
      if not os.path.isfile(fastq):
        raise IOError('Could not find input FASTQ file "'+fastq+'".')
    if args_samples:
      samples = args_samples
    else:
      # Create sample names from fastq filenames.
      samples = []
      for fastq in fastqs1:
        filename = os.path.basename(fastq)
        base = os.path.splitext(filename)[0]
        if base.endswith('_1'):
          samples.append(base[:-2])
        else:
          samples.append(base)
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
    return None, None, None
  return fastqs1, fastqs2, samples


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


def run_pipelines(script_dir, ref, samples, fastqs1, fastqs2, outdirs, pipeline_args, begin=1,
                  end=1000, script=None):
  procs = []
  for sample, fastq1, fastq2, outdir in zip(samples, fastqs1, fastqs2, outdirs):
    if script:
      script_path = os.path.join(outdir, script)
      script_args = '-b '+script_path+' --script-mode a'
    else:
      script_args = ''
    command = ('{script_dir}/pipeline.py {script_dir}/ref-free2.yaml {ref} {fastq1} {fastq2} '
               '{outdir} -B {begin} -E {end} -s {sample} {script_args} {pipeline_args}'
               .format(script_dir=script_dir, ref=ref, fastq1=fastq1, fastq2=fastq2, outdir=outdir,
                       sample=sample, begin=begin, end=end, script_args=script_args,
                       pipeline_args=' '.join(pipeline_args)))
    # Execute command in the background (child process).
    procs.append(subprocess.Popen(command, shell=True))
  # Wait until all samples are done.
  for proc in procs:
    proc.wait()


def parse_report(report_str, report_codes):
  """Parse an asm-curator.py report into a dict.
  Returns a defaultdict with a default value of False, to avoid the mess of KeyErrors."""
  report = collections.defaultdict(bool)
  for line in report_str.splitlines():
    fields = line.strip().split('\t')
    if len(fields) != 3:
      continue
    try:
      key = report_codes[fields[0]]
    except KeyError:
      continue
    value = fields[1]
    if value.lower() == 'true':
      value = True
    elif value.lower() == 'false':
      value = False
    else:
      try:
        value = int(value)
      except ValueError:
        pass
    report[key] = value
  return report


def choose_asm(reports, max_contigs, include_dup):
  """Choose an assembly by prioritizing assembly issues. Some issues are
  dealbreakers, so if all the assemblies have one of these issues, None is
  returned.
  Priority list:
    assembly failure (dealbreaker)
    fragmented assembly (dealbreaker)
    whole-genome duplication (dealbreaker if not include_dup)
    gaps
    non-reference flanks
    number of curated contigs (fewest wins)
    number of original contigs (fewest wins)
  """
  candidates = []
  # Re-score fragmentation according to our limit.
  for report in reports.values():
    if max_contigs and report['contigs-raw'] > max_contigs:
      report['fragmented'] = True
    else:
      report['fragmented'] = False
  # Eliminate samples with dealbreaker issues.
  candidates = reports.keys()
  candidates = asms_without('failure', candidates, reports)
  candidates = asms_without('fragmented', candidates, reports)
  # Whole-genome duplications are dealbreakers unless include_dup.
  if include_dup:
    # If they all have duplication, keep them instead of discarding them.
    candidates = narrow_by('duplication', candidates, reports)
  else:
    # Eliminate any with duplication, even if it's all of them.
    candidates = asms_without('duplication', candidates, reports)
  # Samples without gaps and non-reference flanks are preferred.
  candidates = narrow_by('gaps', candidates, reports)
  candidates = narrow_by('non-ref', candidates, reports)
  # Of the final candidates, choose the one with the fewest curated contigs.
  # If there's a tie, choose the one with the fewest original contigs.
  candidates.sort(key=lambda sample: (reports[sample]['contigs-after'],
                                      reports[sample]['contigs-raw']))
  if candidates:
    return candidates[0]
  else:
    return None


def asms_without(issue, candidates, reports):
  """Return all candidates without the given issue."""
  without = []
  for candidate in candidates:
    if candidate not in reports:
      continue
    if not reports[candidate][issue]:
      without.append(candidate)
  return without


def narrow_by(issue, samples, reports):
  """Use the given issue to narrow the field.
  Eliminate assemblies with the issue, unless all have the issue. Then return
  the original set of assemblies (it's not useful for narrowing the field)."""
  without = asms_without(issue, samples, reports)
  if len(without) > 0:
    return without
  else:
    return samples


def relative_link(src, dst):
  """Make a link to dst at the path src, using a relative path.
  E.g. relative_link('tests/B/file.txt', 'tests/A/file.txt') will make a link at
  'tests/B/file.txt' pointing to '../A/file.txt'."""
  src_dir = os.path.dirname(src)
  rel_path = os.path.relpath(dst, src_dir)
  os.symlink(rel_path, src)


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
