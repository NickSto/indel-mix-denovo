#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess
import distutils.spawn

PLATFORM = 'ILLUMINA'
PICARDIR_DEFAULT = 'bx/src/picard-tools-1.100'
REQUIRED_COMMANDS = ('java', 'bash', 'python', 'bwa', 'samtools', 'spades.py', 'lastz',
                     'naive_variant_caller.py')
REQUIRED_SCRIPTS = ('heteroplasmy/pre-process-mt.sh', 'asm-unifier.py', 'nvc-filter.py',
                    'inspect-reads.py', 'quick-liftover.py')
REQUIRED_PICARDS = ('SamToFastq.jar',)

OPT_DEFAULTS = {'begin':0, 'end':14, 'freq':1, 'cvg':1000, 'strand':1, 'mate':1,
                'margin':600, 'kmers':'21,33,55,77', 'read_length_minimum':40}
DESCRIPTION = """Automate all steps of running the indel discovery pipeline on a single sample.
Note: Set the directory containing the picard jars with the environment variable PICARD_DIR (if it's
different from the default: ~/"""+PICARDIR_DEFAULT+').'
# The intermediate files.
# Their names, without extensions, will be used directly in the command templates below.
# The intention is to make it easy to read the script and correlate filename placeholders in
# commands with the actual ones in the filesystem.
FILENAMES = ('asmdir', 'asm.fa', 'asmlog.log', 'lav.lav', 'bamraw.bam', 'bamfilt.bam',
             'bamdedup.bam', 'nvc.vcf', 'nvcfilt.vcf', 'vars_asm.tsv', 'vars.tsv')


#################### PIPELINE STEPS ####################

# Steps are either a command template that's filled in with parameters and file paths before
# execution by the shell or a call to a function.
STEPS = [
  #TODO: Set k-mers based on read length.
  #TODO: Check if successful (produces a contigs.fasta).
  {'num':1, 'desc':'Assemble with SPAdes', 'type':'command',
   'command':'spades.py --careful -k {kmers} -1 {fastq1} -2 {fastq2} -o {asmdir}'},

  {'num':2, 'desc':'Clean assembly with asm-unifier.py', 'type':'command',
   'command':'python {scriptdir}/asm-unifier.py -n {sample} {ref} {asmdir}/contigs.fasta '
             '-o {asm} -l {asmlog}'},

  {'num':3, 'desc':'Align assembly to reference with LASTZ', 'type':'command',
   'command':'lastz {ref} {asm} > {lav}'},

  {'num':4, 'desc':'Align to assembly with BWA MEM', 'type':'function', 'function':'align',
   'args':['fastq1', 'fastq2', 'asm', 'bamraw', 'sample']},

  # In this second time around, might want to omit -s realign (do NM-edits filtering too).
  #TODO: Check what -c names (name of the assembly sequence) are valid in pre-process-mt.sh.
  {'num':5, 'desc':'Filter alignment with pre-process-mt.sh', 'type':'command',
   'command':'bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {sample} -m {margin} '
             '-M {min_rlen} -s realign -r {asm} {bamraw} {bamfilt}'},

  {'num':6, 'desc':'Remove duplicates with Samtools:', 'type':'function', 'function':'dedup',
   'args':['bamfilt', 'bamdedup']},

  #TODO: Allow setting -q and -m.
  {'num':7, 'desc':'Extract variants with Naive Variant Caller', 'type':'command',
   'command':'naive_variant_caller.py -q 30 -m 20 --ploidy 2 --use_strand '
             '--coverage_dtype uint32 --allow_out_of_bounds_positions --bam {bamdedup} '
             '--index {bamdedup}.bai -r {asm} -o {nvc}'},

  #TODO: Allow for samples with no indels (empty {nvc}).
  #TODO: Postpone application of filters until the end so they can be tweaked. For instance,
  #      we want to preserve information about all indels, even below the MAF threshold.
  {'num':8, 'desc':'Filter variants with nvc-filter.py', 'type':'command',
   'command':'python {scriptdir}/nvc-filter.py -r S -c {cvg} -f {freq} {nvc} > {nvcfilt}'},

  # N.B.: inspect-reads.py takes the output sample name from the BAM read group.
  {'num':9, 'desc':'Filter further by BAM stats with inspect-reads.py', 'type':'command',
   'command':'python {scriptdir}/inspect-reads.py -tl -s {strand} -m {mate} {bamdedup} '
             '-V {nvcfilt} -r {asm} > {vars_asm}'},

  {'num':10, 'desc':'Convert to reference coordinates with quick-liftover.py', 'type':'command',
   'command':'python {scriptdir}/quick-liftover.py {lav} {vars_asm} > {vars}'},
]


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

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
  parser.add_argument('-n', '--simulate', action='store_true',
    help='Only simulate execution. Print commands, but do not execute them.')
  parser.add_argument('-b', '--to-script', metavar='path/to/script.sh',
    help='Instead of executing commands, write them to a bash script with this name.')
  parser.add_argument('-B', '--begin', metavar='step', type=int,
    help='Start at this step. The input files are the normal intermediate files generated if the '
         'full pipeline were run with the same arguments. WARNING: Any existing intermediate files '
         'in the output directory from later steps will be overwritten.')
  parser.add_argument('-E', '--end', metavar='step', type=int,
    help='Stop after this many pipeline steps.')
  param = parser.add_argument_group('Analysis Parameters')
  param.add_argument('-l', '--read-length', dest='rlen', required=True, type=int,
    help='Read length. Note: all reads don\'t have to be this length (variable length reads are '
         'accepted). Default: "%(default)s"')
  param.add_argument('-L', '--read-length-minimum', metavar='PCT', type=float,
    help='Read length threshold. Give the minimum percent of the read that must be present to '
         'pass. Give in percent, not decimal ("10" for 10%%, not "0.1"). Default: "%(default)s"')
  param.add_argument('-f', '--freq-thres', dest='freq', type=float,
    help='Minor allele frequency threshold for indel calling. Give in percent, not decimal. Used '
         'in step 12 (nvc-filter.py). Default: "%(default)s"')
  param.add_argument('-c', '--cvg-thres', dest='cvg', type=int,
    help='Read depth of coverage threshold for indel calling. If the read depth at the indel is '
         'below this value, it will not be reported. Used in step 12 (nvc-filter.py). Default: '
         '"%(default)s"')
  param.add_argument('-S', '--strand-bias', dest='strand', type=float,
    help='Strand bias threshold. Used in step 13 (inspect-reads.py). Default: "%(default)s"')
  param.add_argument('-M', '--mate-bias', dest='mate', type=float,
    help='Mate bias threshold. Used in step 13 (inspect-reads.py). Default: "%(default)s"')
  param.add_argument('-k', '--k-mers', dest='kmers',
    help='K-mers for assembly. Comma-delimited list of ascending integers. Will be passed directly '
         'to spades.py. Default: "%(default)s"')
  param.add_argument('-m', '--margin', type=int,
    help='Size of the regions at either end of the reference where chimeric reads are allowed. '
         'Used for circular chromosomes where reads spanning the start coordinate appear as '
         'chimeric but aren\'t. Give a size in nucleotides. Both margins, at the start and end, '
         'will be this size. Set to 0 for no margins. Default: "%(default)s"')
  args = parser.parse_args(argv[1:])

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

  # Create a Runner and set it to execute the commands immediately, just print them to stdout
  # (simulate), or output them to a script file.
  runner = Runner()
  runner.simulate = args.simulate
  if args.to_script:
    runner.to_script(args.to_script)

  # Create paths to files and directories.
  #TODO: Use a tmp directory for intermediate files
  params = {}
  for filename in FILENAMES:
    base = os.path.splitext(filename)[0]
    assert base not in params, '{} in params. value: {}'.format(base, params[base])
    params[base] = os.path.join(args.outdir, filename)
  params['scriptdir'] = os.path.dirname(os.path.realpath(sys.argv[0]))
  if 'PICARD_DIR' in os.environ:
    params['picardir'] = os.environ['PICARD_DIR']
  else:
    params['picardir'] = os.path.join(os.path.expanduser('~'), PICARDIR_DEFAULT)
  # Add arguments from command line.
  for arg in dir(args):
    if arg.startswith('_'):
      continue
    assert arg not in params, '{} in params. value: {}'.format(arg, params[arg])
    params[arg] = getattr(args, arg)
  # Compute some special arguments.
  if not params.get('rlen'):
    #TODO: Detect read length of FASTA files.
    raise Exception('Need to provide --read-length.')
  params['min_rlen'] = int(round(params['rlen'] * params['read_length_minimum'] / 100))

  # Check for required commands.
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      raise Exception('Required command "'+command+'" not found.')
  # Check for required scripts.
  for script in REQUIRED_SCRIPTS:
    scriptpath = os.path.join(params['scriptdir'], script)
    if not os.path.isfile(scriptpath):
      raise Exception('Required script "'+scriptpath+'" not found.')
  # Check for Picard jars.
  for jar in REQUIRED_PICARDS:
    jarpath = os.path.join(params['picardir'], jar)
    if not os.path.isfile(jarpath):
      raise Exception('Required Picard jar "'+jarpath+'" not found.')


  #################### EXECUTE STEPS ####################

  # Get a dict of the global variables so we can look up functions by name.
  globals_dict = globals()
  for step in STEPS:
    if args.begin <= step['num'] and args.end >= step['num']:
      print 'Step {num}: {desc}.'.format(**step)
      if step['type'] == 'command':
        # Fill in the placeholders in the command and have the Runner "execute" it.
        runner.run(step['command'].format(**params))
      elif step['type'] == 'function':
        # Look up the requested function and its arguments, then call it.
        function = globals_dict[step['function']]
        arguments = [params[arg] for arg in step['args']]
        function(runner, *arguments)
    elif args.end > step['num']:
      print 'Skipping step {num}: {desc}.'.format(**step)

  if args.end >= STEPS[-1]['num']:
    print 'Finished. Final variants are in {vars}.'.format(**params)


def align(runner, fastq1, fastq2, ref, outbam, sample):
  """Map the reads in "fastq1" and "fastq2" to reference "ref", output to "outbam". "sample" will be
  used to label read groups."""
  if runner is None:
    runner = Runner()
  (base, ext) = os.path.splitext(outbam)
  if ext != '.bam':
    base = base+ext
  rg_line = "'@RG\\tID:{}\\tSM:{}\\tPL:{}'".format(sample, sample, PLATFORM)
  params = {'rg':rg_line, 'ref':ref, 'fq1':fastq1, 'fq2':fastq2, 'outbam':outbam, 'base':base}
  # Index reference, if needed.
  for ext in ('.amb', '.ann', '.bwt', '.sa', '.pac'):
    if not os.path.isfile(ref+ext):
      #TODO: determine algorithm based on file size ("is" if size is less than 2000000000 bytes).
      algorithm = 'bwtsw'
      runner.run('bwa index -a {algo} {ref}'.format(algo=algorithm, ref=ref))
      break
  runner.run('bwa mem -M -t 16 -R {rg} {ref} {fq1} {fq2} > {base}.sam'.format(**params))
  runner.run('samtools view -Sb {base}.sam > {base}.tmp.bam'.format(**params))
  runner.run('samtools sort {base}.tmp.bam {base}'.format(**params))
  #TODO: fixmate?
  runner.run('samtools index {base}.bam'.format(**params))
  runner.run('rm {base}.sam {base}.tmp.bam'.format(**params))
  return outbam


def dedup(runner, inbam, outbam):
  (base, ext) = os.path.splitext(outbam)
  params = {'inbam':inbam, 'base':base}
  runner.run('samtools view -b -F 1024 {inbam} > {base}.tmp.bam'.format(**params))
  runner.run('samtools sort {base}.tmp.bam {base}'.format(**params))
  runner.run('samtools index {base}.bam'.format(**params))
  runner.run('rm {base}.tmp.bam'.format(**params))


class Runner(object):
  """An object that will execute commands according to previously arranged settings."""
  def __init__(self):
    self.quiet = False
    self.simulate = False
    self.prepend = '+ '
    self._output = sys.stdout
  def run(self, command, ignore_err=False):
    if not self.quiet:
      self._output.write(self.prepend+command+'\n')
    if not self.simulate:
      try:
        subprocess.check_call(command, shell=True)
      except subprocess.CalledProcessError as cpe:
        sys.stderr.write("Command '{}' returned non-zero exit status {}\n"
                         .format(cpe.cmd, cpe.returncode))
        if not ignore_err:
          if __name__ == '__main__':
            sys.exit(cpe.returncode)
          else:
            raise
  def to_script(self, path):
    """Print the commands to a ready-to-run bash script instead of executing them."""
    self.quiet = False
    self.simulate = True
    self.prepend = ''
    self._output = open(path, 'w')
  #TODO: close filehandle on exit or when done


if __name__ == '__main__':
  sys.exit(main(sys.argv))
