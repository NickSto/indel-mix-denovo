#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess
import distutils.spawn
import yaml

PLATFORM = 'ILLUMINA'
PICARD_DIR_DEFAULT = '~/bx/src/picard-tools-1.100'

DESCRIPTION = """Automate all steps of running the indel discovery pipeline on a single sample."""

# About "filenames" in YAML files:
# You can use these names, without the extensions, directly in command templates.
# The names, with extensions, will be substituted into the executed command.
# The intention is to make it easy to read the script and correlate filename placeholders in
# commands with the actual ones in the filesystem.


#################### PIPELINE STEPS ####################


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('yaml',
    help='A YAML file with the pipeline commands.')
  parser.add_argument('ref', metavar='reference.fa',
    help='The reference genome.')
  parser.add_argument('fastq1', metavar='reads_1.fq',
    help='Input reads, mate 1.')
  parser.add_argument('fastq2', metavar='reads_2.fq',
    help='Input reads, mate 2.')
  parser.add_argument('outdir', metavar='output/directory/path',
    help='Destination directory to place the output.')
  parser.add_argument('-s', '--sample', metavar='id', required=True,
    help='The sample id. Required.')
  parser.add_argument('--refname', metavar='chrName',
    help='The id of the reference sequence.')
  parser.add_argument('--refname2', metavar='chrName',
    help='Id of the sequence to filter for, after assembly (if different from the --sample). '
         'Normally only needed when starting the pipeline at a middle step after changing the '
         'assembly file. If not given, will default to the --sample.')
  parser.add_argument('-R', '--filter-ref', metavar='filter-ref.fa',
    help='A reference consisting of the reference genome plus any other mapping targets you want '
         'to include for filtering purposes. Will filter the reads in the first step by mapping '
         'to this reference, then selecting only reads which map to the sequence identified by '
         '--refname. Optional. Will use the main reference genome by default.')
  parser.add_argument('-P', '--picard-dir',
    help='The directory containing Picard jar files. Default: the value of $PICARD_DIR, or '+
          PICARD_DIR_DEFAULT+' if unset. Using this option sets $PICARD_DIR for this script and '
          'its children.')
  parser.add_argument('-n', '--simulate', dest='execute', action='store_false', default=True,
    help='Only simulate execution. Print commands, but do not execute them.')
  parser.add_argument('-b', '--to-script', metavar='path/to/script.sh',
    help='Write the commands to this bash script.')
  parser.add_argument('--script-mode', choices=('w', 'a'), default='w',
    help='The file mode to write to the --to-script. "w" for write, "a" for append. Default: '
         '%(default)s')
  parser.add_argument('-B', '--begin', metavar='step', type=int, default=0,
    help='Start at this step. The input files are the normal intermediate files generated if the '
         'full pipeline were run with the same arguments. WARNING: Any existing intermediate files '
         'in the output directory from later steps will be overwritten.')
  parser.add_argument('-E', '--end', metavar='step', type=int, default=14,
    help='Stop after this many pipeline steps.')
  parser.add_argument('-t', '--threads', type=int,
    help='Number of threads to use for tools like bwa and SPAdes. Default: %(default)s')
  param = parser.add_argument_group('Analysis Parameters')
  param.add_argument('-l', '--read-length', dest='rlen', required=True, type=int,
    help='Read length. Note: all reads don\'t have to be this length (variable length reads are '
         'accepted). Default: "%(default)s"')
  param.add_argument('-L', '--read-length-minimum', metavar='PCT', type=float, default=40,
    help='Read length threshold. Give the minimum percent of the read that must be present to '
         'pass. Give in percent, not decimal ("10" for 10%%, not "0.1"). Default: "%(default)s"')
  param.add_argument('-f', '--freq-thres', dest='freq', type=float, default=1,
    help='Minor allele frequency threshold for indel calling. Give in percent, not decimal. Used '
         'in step 12 (nvc-filter.py). Default: "%(default)s"')
  param.add_argument('-c', '--cvg-thres', dest='cvg', type=int, default=1000,
    help='Read depth of coverage threshold for indel calling. If the read depth at the indel is '
         'below this value, it will not be reported. Used in step 12 (nvc-filter.py). Default: '
         '"%(default)s"')
  param.add_argument('-S', '--strand-bias', dest='strand', type=float, default=1,
    help='Strand bias threshold. Used in step 13 (inspect-reads.py). Default: "%(default)s"')
  param.add_argument('-M', '--mate-bias', dest='mate', type=float, default=1,
    help='Mate bias threshold. Used in step 13 (inspect-reads.py). Default: "%(default)s"')
  param.add_argument('-e', '--excluded', metavar="302-310,chrM:16183-16192",
    help='Regions to exclude. Variants in these regions (including the ends) will be filtered out. '
         'Chromosome names are optional; if not given, the region will be excluded from every '
         'chromosome. N.B.: chromosome names can\'t contain a dash.')
  param.add_argument('-k', '--k-mers', dest='kmers', default='21,33,55,77',
    help='K-mers for assembly. Comma-delimited list of ascending integers. Will be passed directly '
         'to spades.py. Default: "%(default)s"')
  param.add_argument('-m', '--margin', type=int, default=600,
    help='Size of the regions at either end of the reference where chimeric reads are allowed. '
         'Used for circular chromosomes where reads spanning the start coordinate appear as '
         'chimeric but aren\'t. Give a size in nucleotides. Both margins, at the start and end, '
         'will be this size. Set to 0 for no margins. Default: "%(default)s"')
  args = parser.parse_args(argv[1:])

  with open(args.yaml) as yaml_file:
    script = yaml.safe_load(yaml_file)

  # Check output directory.
  if args.begin > 1:
    # If beginning later in the pipeline, the outdir must exist and contain needed files.
    if not os.path.isdir(args.outdir):
      raise Exception('Output directory "'+args.outdir+'" missing. Must be present and contain '
                      'intermediate files if using --begin.')
  else:
    # Normally, output directory must either not exist or be empty.
    if os.path.exists(args.outdir):
      if os.path.isdir(args.outdir):
        if os.listdir(args.outdir):
          raise Exception('Output directory "'+args.outdir+'" exists but is not empty.')
      else:
        raise Exception('Output directory "'+args.outdir+'" exists but is not a directory.')
    else:
      # Make the output directory if it doesn't exist.
      os.makedirs(args.outdir)

  # Create a Runner, with the user-requested settings.
  runner = Runner(execute=args.execute)
  if args.to_script:
    runner.script = open(args.to_script, args.script_mode)
    if args.script_mode == 'w':
      runner.script.write('#!/usr/bin/env bash\n')

  # Gather all parameters: file paths from filenames in YAML, command line arguments, etc.
  # File paths.
  #TODO: Use a tmp directory for intermediate files
  params = {}
  for filename in script['filenames']:
    base = os.path.splitext(filename)[0]
    assert base not in params, '{} in params. value: {}'.format(base, params[base])
    params[base] = os.path.join(args.outdir, filename)
  params['scriptdir'] = os.path.dirname(os.path.realpath(sys.argv[0]))
  if args.picard_dir:
    params['picardir'] = args.picard_dir
    os.environ['PICARD_DIR'] = args.picard_dir
  elif 'PICARD_DIR' in os.environ:
    params['picardir'] = os.environ['PICARD_DIR']
  else:
    params['picardir'] = os.path.expanduser(PICARD_DIR_DEFAULT)
  # Command line arguments.
  for arg in dir(args):
    if arg.startswith('_'):
      continue
    # Command line argument names can't conflict with filenames or other existing parameters.
    assert arg not in params, '{} in params. value: {}'.format(arg, params[arg])
    params[arg] = getattr(args, arg)
  # Compute some special parameters.
  if not args.filter_ref:
    params['filter_ref'] = args.ref
  if not args.refname2:
    params['refname2'] = args.sample
  if not params.get('rlen'):
    #TODO: Detect read length of FASTQ files.
    raise Exception('Need to provide --read-length.')
  params['min_rlen'] = int(round(params['rlen'] * params['read_length_minimum'] / 100))

  # Check for required commands.
  requirements = script.get('requirements', {})
  for command in requirements.get('commands', ()):
    if not distutils.spawn.find_executable(command):
      raise Exception('Required command "'+command+'" not found.')
  # Check for required scripts.
  for script_name in requirements.get('scripts', ()):
    scriptpath = os.path.join(params['scriptdir'], script_name)
    if not os.path.isfile(scriptpath):
      raise Exception('Required script "'+scriptpath+'" not found.')
  # Check for Picard jars.
  for jar in requirements.get('picards', ()):
    jarpath = os.path.join(params['picardir'], jar)
    if not os.path.isfile(jarpath):
      raise Exception('Required Picard jar "'+jarpath+'" not found.')


  #################### EXECUTE STEPS ####################

  # Steps are either a command template that's filled in with parameters and file paths before
  # execution by the shell or a call to a function.

  # Get a dict of the global variables so we can look up functions by name.
  globals_dict = globals()
  for step in script['steps']:
    if args.begin <= step['num'] and args.end >= step['num']:
      print 'Step {num}: {desc}.'.format(**step)
      if runner.script:
        runner.script.write('# Step {num}: {desc}.\n'.format(**step))
      if 'command' in step:
        kwargs = {}
        if 'ignore_err' in step:
          kwargs['ignore_err'] = step['ignore_err']
        # Fill in the placeholders in the command and have the Runner "execute" it.
        runner.run(step['command'].format(**params), **kwargs)
      elif 'function' in step:
        # Look up the requested function and its arguments, then call it.
        function = globals_dict[step['function']]
        arguments = [params[arg] for arg in step['args']]
        function(runner.run, *arguments)
    elif args.end > step['num']:
      print 'Skipping step {num}: {desc}.'.format(**step)

  last_step = script['steps'][-1]
  if args.end >= last_step['num']:
    print 'Finished. Final variants are in {output}.'.format(**script)


def align(runfxn, fastq1, fastq2, ref, outbam, sample, threads):
  """Map the reads in "fastq1" and "fastq2" to reference "ref", output to "outbam". "sample" will be
  used to label read groups."""
  if runfxn is None:
    runfxn = Runner().run
  (base, ext) = os.path.splitext(outbam)
  if ext != '.bam':
    base = base+ext
  rg_line = "'@RG\\tID:{}\\tSM:{}\\tPL:{}'".format(sample, sample, PLATFORM)
  params = {
    'rg':rg_line, 'ref':ref, 'fq1':fastq1, 'fq2':fastq2, 'outbam':outbam, 'base':base,
    'threads':threads,
  }
  # Index reference, if needed.
  for ext in ('.amb', '.ann', '.bwt', '.sa', '.pac'):
    if not os.path.isfile(ref+ext):
      #TODO: determine algorithm based on file size ("is" if size is less than 2000000000 bytes).
      algorithm = 'bwtsw'
      runfxn('bwa index -a {algo} {ref}'.format(algo=algorithm, ref=ref))
      break
  runfxn('bwa mem -M -t {threads} -R {rg} {ref} {fq1} {fq2} > {base}.sam'.format(**params))
  runfxn('samtools view -Sb {base}.sam > {base}.tmp.bam'.format(**params))
  runfxn('samtools sort {base}.tmp.bam {base}'.format(**params))
  #TODO: fixmate?
  runfxn('samtools index {base}.bam'.format(**params))
  runfxn('rm {base}.sam {base}.tmp.bam'.format(**params))
  return outbam


def dedup(runfxn, inbam, outbam):
  if runfxn is None:
    runfxn = Runner().run
  (base, ext) = os.path.splitext(outbam)
  params = {'inbam':inbam, 'base':base}
  runfxn('samtools view -b -F 1024 {inbam} > {base}.tmp.bam'.format(**params))
  runfxn('samtools sort {base}.tmp.bam {base}'.format(**params))
  runfxn('samtools index {base}.bam'.format(**params))
  runfxn('rm {base}.tmp.bam'.format(**params))


class Runner(object):
  """An object that will execute commands according to previously arranged settings."""
  def __init__(self, execute=True, echo=True, prepend='+ ', stream=sys.stdout, script=None):
    self.execute = execute
    self.echo = echo
    self.prepend = prepend
    self.stream = stream
    self.script = script
  def run(self, command, ignore_err=False):
    command_stripped = command.strip()
    if self.echo:
      self.stream.write(self.prepend+command_stripped+'\n')
    if self.script:
      self.script.write(command_stripped+'\n')
    if self.execute:
      try:
        subprocess.check_call(command_stripped, shell=True)
      except subprocess.CalledProcessError as cpe:
        sys.stderr.write("Command '{}' returned non-zero exit status {}\n"
                         .format(cpe.cmd, cpe.returncode))
        if not ignore_err:
          if __name__ == '__main__':
            sys.exit(cpe.returncode)
          else:
            raise


if __name__ == '__main__':
  sys.exit(main(sys.argv))
