#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
import time
import copy
import shutil
import argparse
import subprocess
import collections
import multiprocessing

OPT_DEFAULTS = {'max_contigs':500}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
FASTQ_NAME_REGEX = r'(.+)_[12]\.(fastq|fq)(\.gz)?$'

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
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('ref', metavar='reference.fa',
    help='The reference genome.')
  parser.add_argument('outdir', metavar='output/directory/path',
    help='Destination directory to place the output.')
  parser.add_argument('-r', '--refname', metavar='chrName', required=True,
    help='Id of reference sequence. Reads will be filtered for those that map to this sequence.')
  parser.add_argument('-l', '--read-length', required=True, type=int,
    help='Read length. Default: "%(default)s"')
  parser.add_argument('-f', '--fastq-dir', required=True,
    help='Directory containing the fastq files. The filenames should be in the format '
         '{sampleid}_1.fastq.gz.')
  parser.add_argument('-t', '--families-table', required=True,
    help='A file containing the list of families. Format: A tab-delimited file with one line per '
         'each family. The first column is the id of the family, and the following columns are the '
         'ids of the samples in that family.')
  parser.add_argument('-F', '--family-id')
  parser.add_argument('-D', '--include-dup', action='store_true',
    help='Don\'t fail if all assemblies in the family show whole-genome duplications.')
  parser.add_argument('-C', '--max-contigs',
    help='The maximum allowed number of contigs per assembly (approximated by the number of LASTZ '
         'hits to the reference). Set to 0 to allow any number. Default: %(default)s')
  parser.add_argument('-N', '--simulate', action='store_true',
    help='Only simulate execution. Print commands, but do not execute them.')
  parser.add_argument('-b', '--to-script', metavar='path/to/script.sh',
    help='Instead of executing commands, write them to a bash script with this name.')
  parser.add_argument('pipeargs', metavar='...', nargs=argparse.REMAINDER,
    help='The remaining arguments will be passed directly to pipeline.py. Make sure to put them '
         'at the end, after the arguments for this script.')
  parser.add_argument('-p', '--pre', metavar='cmd -arg1',
    help='A command to prepend to all shell commands. For instance, "-p srun" to run the commands '
         'on the SLURM job manager. You can give arguments to the command too; just put the whole '
         'thing in a quoted string, like so: "-p \'srun -C grp1\'". This will transform a command '
         'line like "$ python pipeline.py -s 1" into "$ srun -C grp1 python pipeline.py -s 1".')

  args = parser.parse_args(argv[1:])

  pre = []
  if args.pre:
    pre = args.pre.split()

  script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

  fastq_pairs = get_fastqs(args.fastq_dir)

  families = read_families(args.families_table)
  for family_id, sample_ids in families.items():
    print 'family '+family_id
    # Run first half of pipeline (up through assembly).
    processes = []
    outdirs = {}
    for sample_id in sample_ids:
      if sample_id in fastq_pairs:
        fastq_pair = fastq_pairs[sample_id]
        if len(fastq_pair) != 2:
          raise PipefamError('Wrong number of fastq\'s associated with sample "{}": {}'
                             .format(sample_id, ', '.join(fastq_pair)))
      else:
        raise PipefamError('No fastq\'s found for sample ""'.format(sample_id))
      # Make output directory.
      outdirs[sample_id] = os.path.join(args.outdir, sample_id)
      makedir_and_check(outdirs[sample_id])
      # Build pipeline.py command.
      fastq1 = os.path.join(args.fastq_dir, fastq_pair[0])
      fastq2 = os.path.join(args.fastq_dir, fastq_pair[1])
      command = copy.copy(pre)
      command.extend(['python', os.path.join(script_dir, 'pipeline.py'), '-s', sample_id,
                      '-r', args.refname, '-l', str(args.read_length), '-E', '7'])
      command.extend(args.pipeargs)
      command.extend([args.ref, fastq1, fastq2, outdirs[sample_id]])
      # command.extend(['-b', '/dev/null'])
      print '+ $ '+' '.join(command)
      if not args.simulate:
        process = multiprocessing.Process(target=subprocess.call, args=(command,))
        process.start()
        processes.append(process)
    while not processes_done(processes):
      time.sleep(1)
    print "Family {} assembled.".format(family_id)
    # Check assemblies and choose the best one.
    reports = {}
    for sample_id in sample_ids:
      asm_raw = os.path.join(outdirs[sample_id], 'asmdir', 'contigs.fasta')
      if not os.path.isfile(asm_raw):
        continue
      # Align raw assembly to reference.
      lav_path = os.path.join(outdirs[sample_id], 'lav_raw.lav')
      command = copy.copy(pre)
      command.extend(['lastz', args.ref, asm_raw])
      with open(lav_path, 'w') as lav:
        subprocess.call(command, stdout=lav)
      if not os.path.isfile(lav_path) or os.path.getsize(lav_path) == 0:
        continue
      # Use asm-curator.py to generate statistics.
      command = copy.copy(pre)
      command.extend(['python', os.path.join(script_dir, 'asm-curator.py'), '-l', lav_path])
      print '+ $ '+' '.join(command)
      if not args.simulate:
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        reports[sample_id] = parse_report(subprocess.communicate()[0])
    # Compare reports for all samples and choose the best assembly.
    if not reports:
      sys.stderr.write('No successful assemblies found for family "{}".\n'.format(family_id))
      continue
    best_sample = choose_asm(reports, args.max_contigs, args.include_dup)
    if best_sample is None:
      sys.stderr.write('No high-quality assembly chosen for family "{}".\n'.format(family_id))
      continue
    best_asm = os.path.join(outdirs[best_sample], 'asm.fa')
    if not os.path.isfile(best_asm):
      continue
    # Finish pipeline using one, best assembly for all samples.
    processes = []
    for sample_id in sample_ids:
      asm = os.path.join(outdirs[sample_id], 'asm.fa')
      #TODO: Catch shutil.Error and OSError (or IOError?)
      shutil.copy2(best_asm, asm)
      # Build pipeline.py command.
      fastq_pair = fastq_pairs[sample_id]
      fastq1 = os.path.join(args.fastq_dir, fastq_pair[0])
      fastq2 = os.path.join(args.fastq_dir, fastq_pair[1])
      command = copy.copy(pre)
      command.extend(['python', os.path.join(script_dir, 'pipeline.py'), '-s', sample_id,
                      '-r', args.refname, '-l', str(args.read_length), '-B', '8'])
      print '+ $ '+' '.join(command)
      if not args.simulate:
        process = multiprocessing.Process(target=subprocess.call, args=(command,))
        process.start()
        processes.append(process)
    while not processes_done(processes):
      time.sleep(1)
    print "Family {} done!".format(family_id)
    #TODO: Consider bringing back group-filter.py method of allowing indels above bias thresholds
    #      if there are others in the family that pass it.


def read_families(families_filepath):
  """Read the families tsv file and return a dict mapping family ID's to lists of sample ID's in
  that family."""
  families = {}
  with open(families_filepath) as families_file:
    for line in families_file:
      fields = line.rstrip('\r\n').split('\t')
      if len(fields) < 2:
        continue
      family_id = None
      sample_ids = []
      for field in fields:
        if family_id is None:
          family_id = field
        else:
          sample_ids.append(field)
      families[family_id] = sample_ids
  return families


def get_fastqs(fastq_dir):
  """Read the files in a directory and return a mapping of sample ID's to fastq filenames.
  The fastqs are found by matching FASTQ_NAME_REGEX and sample_id's are extracted that way.
  Returns a dict mapping each sample_id to a list of the filenames of fastq's that match that
  sample_id."""
  fastqs = {}
  # Read the files, build dict of fastqs.
  for filename in os.listdir(fastq_dir):
    match = re.search(FASTQ_NAME_REGEX, filename)
    if not match:
      return
    sample_id = match.group(1)
    fastq_pair = fastqs.get(sample_id, [])
    fastq_pair.append(filename)
    fastqs[sample_id] = fastq_pair
  # Check and sort the fastq filename lists.
  for sample_id, fastq_pair in fastqs.items():
    if len(fastq_pair) != 2:
      raise PipefamError('More or less than 2 files match the same sample ID "{}": {}'
                         .format(sample_id, ', '.join(fastq_pair)))
    fastq_pair.sort()
  return fastqs


def processes_done(processes):
  """Return True if all processes are finished, False otherwise."""
  for process in processes:
    if process.is_alive():
      return False
  return True


def makedir_and_check(dirpath):
  """Make the directory if it does not exist already.
  Raise exception if the path already exists and a) is not a directory or b) is not empty."""
  if os.path.exists(dirpath) and not os.path.isdir(dirpath):
    raise PipefamError('Output directory path '+dirpath+' exists but is not a directory.')
  elif os.path.isdir(dirpath) and os.listdir(dirpath):
    raise PipefamError('Output directory '+dirpath+' is not empty.')
  elif not os.path.exists(dirpath):
    os.makedirs(dirpath)


def parse_report(report_str):
  """Parse an asm-curator.py report into a dict.
  Returns a defaultdict with a default value of False, to avoid the mess of KeyErrors."""
  report = collections.defaultdict(bool)
  for line in report_str.splitlines():
    fields = line.strip().split('\t')
    if len(fields) != 3:
      continue
    try:
      key = REPORT_CODES[fields[0]]
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


class PipefamError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
