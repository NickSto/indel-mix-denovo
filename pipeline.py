#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess

PLATFORM = 'ILLUMINA'

OPT_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
FILENAMES = ('bam1raw.bam', 'bam1filt.bam')


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('ref', metavar='reference.fa',
    help='Reference genome, including any other mapping targets you want to include for filtering '
         'purposes.')
  parser.add_argument('fastq1', metavar='reads_1.fq',
    help='Input reads, mate 1.')
  parser.add_argument('fastq2', metavar='reads_2.fq',
    help='Input reads, mate 2.')
  parser.add_argument('outdir', metavar='output/directory/path',
    help='Destination directory to place the output.')
  parser.add_argument('-s', '--sample', required=True,
    help='The sample id.')
  parser.add_argument('-R', '--refname', required=True,
    help='Id of reference sequence. Reads will be filtered for those that map to this sequence.')
  parser.add_argument('-N', '--simulate', action='store_true',
    help='Only simulate execution. Print commands, but do not execute them.')
  parser.add_argument('-b', '--to-script',
    help='Instead of executing commands, write them to a bash script with this name.')

  args = parser.parse_args(argv[1:])

  if not os.path.isdir(args.outdir):
    raise Exception('Output directory "'+args.outdir+'" is not accessible or not a directory.')

  runner = Runner()
  runner.simulate = args.simulate
  if args.to_script:
    runner.to_script(args.to_script)

  # Create paths to files and directories
  #TODO: Use a tmp directory for intermediate files
  cmd_args = {}
  for filename in FILENAMES:
    base = os.path.splitext(filename)[0]
    cmd_args[base] = os.path.join(args.outdir, filename)
  cmd_args['ref'] = args.ref
  cmd_args['scriptdir'] = os.path.relpath(os.path.dirname(os.path.realpath(sys.argv[0])))
  # Set general arguments for commands
  cmd_args['sample'] = args.sample
  cmd_args['refname'] = args.refname

  # Map reads to reference.
  align(args.fastq1, args.fastq2, args.ref, cmd_args['bam1raw'], args.sample, runner)

  # Filter alignment
  #TODO: only use -s realign when necessary
  #TODO: allow setting margin
  # Example command:
  # $ pre-process-mt.sh -r $root/asm/clean/G3825.1a.fa -c G3825.1a -B '0 18834' -s realign \
  #     $root/aln/raw/G3825.1a.bam $root/aln/filt/G3825.1a.bam $root/aln/tmp/G3825.1a &
  cmd_args['margin'] = 600
  runner.run('bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {refname} -m {margin} -s realign '
             '-r {ref} {bam1raw} {bam1filt}'.format(**cmd_args))

  # Remove duplicates
  #   samtools

  # Extract reads
  #   Picard SamToFastq

  # Assemble
  #   SPAdes

  # Clean assembly
  #   asm-unifier.py

  # Align assembly to reference
  #   LASTZ

  # Align to assembly
  #   BWA-MEM

  # Filter alignment
  #   pre-process-mt.sh

  # Remove duplicates
  #   samtools

  # Naive Variant Caller

  # nvc-filter.py

  # inspect-reads.py

  # quick-liftover.py

  # Final indel filtering
  #   awk


def align(fastq1, fastq2, ref, outbam, sample, runner):
  """Map the reads in "fastq1" and "fastq2" to reference "ref", output to "outbam". "sample" will be
  used to label read groups."""
  if runner is None:
    runner = Runner()
  (base, ext) = os.path.splitext(outbam)
  if ext != '.bam':
    base = base+ext
  rg_line = "'@RG\\tID:{}\\tSM:{}\\tPL:{}'".format(sample, sample, PLATFORM)
  paths = {'rg':rg_line, 'ref':ref, 'fq1':fastq1, 'fq2':fastq2, 'bam':outbam, 'base':base,
           'sam':base+'.sam', 'bamtmp':base+'.tmp.bam'}
  # Index reference, if needed.
  for ext in ('.amb', '.ann', '.bwt', '.sa', '.pac'):
    if not os.path.isfile(ref+ext):
      #TODO: determine algorithm based on file size ("is" if size is less than 2000000000 bytes).
      algorithm = 'bwtsw'
      runner.run('bwa index -a {algo} {ref}'.format(algo=algorithm, ref=ref))
      break
  runner.run('bwa mem -M -t 16 -R {rg} {ref} {fq1} {fq2} > {sam}'.format(**paths))
  runner.run('samtools view -Sb {sam} > {bamtmp}'.format(**paths))
  runner.run('samtools sort {bamtmp} {base}'.format(**paths))
  runner.run('samtools index {bam}'.format(**paths))
  runner.run('rm {sam} {bamtmp}'.format(**paths))
  return outbam


class Runner(object):
  """An object that will execute commands according to previously arranged settings."""
  def __init__(self):
    self.quiet = False
    self.simulate = False
    self.prepend = '+ '
    self._output = sys.stdout
  def run(self, command):
    if not self.quiet:
      self._output.write(self.prepend+command+'\n')
    if not self.simulate:
      subprocess.check_call(command, shell=True)
  def to_script(self, path):
    """Print the commands to a ready-to-run bash script instead of executing them."""
    self.quiet = False
    self.simulate = True
    self.prepend = ''
    self._output = open(path, 'w')
  #TODO: close filehandle on exit or when done


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
