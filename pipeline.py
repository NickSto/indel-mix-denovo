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

  args = parser.parse_args(argv[1:])

  if not os.path.isdir(args.outdir):
    raise Exception('Output directory "'+args.outdir+'" is not accessible or not a directory.')

  runner = Runner()
  runner.simulate = args.simulate

  # Map reads to reference.
  bam = os.path.join(args.outdir, args.sample+'.bam')
  align(args.fastq1, args.fastq2, args.ref, bam, args.sample, runner=runner)


class Runner(object):
  def __init__(self):
    self.quiet = False
    self.simulate = False
    self.prepend = ''

  def run(self, command):
    if not self.quiet:
      print '+ '+command
    if not self.simulate:
      subprocess.check_call(self.prepend+command)


def align(fastq1, fastq2, ref, outbam, sample, runner=Runner()):
  """Map the reads in "fastq1" and "fastq2" to reference "ref", output to "outbam". "sample" will be
  used to label read groups."""
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


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
