#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess

PLATFORM = 'ILLUMINA'
PICARDIR = 'src/picard-tools-1.100'

OPT_DEFAULTS = {'begin':0, 'end':15}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
FILENAMES = ('bam1raw.bam', 'bam1filt.bam', 'bam1dedup.bam', 'cleanfq1.fq', 'cleanfq2.fq',
             'asmdir', 'asm.fa', 'asmlog.log', 'lav.lav', 'bam2raw.bam', 'bam2filt.bam',
             'bam2dedup.bam', 'nvc.vcf', 'nvcfilt.vcf', 'vars1asm.tsv', 'vars1.tsv')


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
  parser.add_argument('-B', '--begin', type=int,
    help='Start at this step. The input files are the normal intermediate files generated if the '
         'full pipeline were run with the same arguments. WARNING: Any existing intermediate files '
         'in the output directory from later steps will be overwritten.')
  parser.add_argument('-E', '--end', type=int,
    help='Stop after this many pipeline steps.')

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

  runner = Runner()
  runner.simulate = args.simulate
  if args.to_script:
    runner.to_script(args.to_script)

  #TODO: Check for required commands.

  # Create paths to files and directories
  #TODO: Use a tmp directory for intermediate files
  paths = {}
  for filename in FILENAMES:
    base = os.path.splitext(filename)[0]
    paths[base] = os.path.join(args.outdir, filename)
  paths['ref'] = args.ref
  paths['scriptdir'] = os.path.relpath(os.path.dirname(os.path.realpath(sys.argv[0])))
  paths['picardir'] = os.path.join(os.path.expanduser('~'), PICARDIR)
  # Set general arguments for commands
  paths['sample'] = args.sample
  paths['refname'] = args.refname

  # Map reads to reference.
  #   BWA MEM
  if args.begin <= 1 and args.end >= 1:
    align(args.fastq1, args.fastq2, args.ref, paths['bam1raw'], args.sample, runner)
  else:
    print 'Skipping step 1.'

  # Filter alignment
  #   pre-process-mt.sh
  #TODO: only use -s realign when necessary
  #TODO: allow setting margin
  # Example command: $ pre-process-mt.sh -r $root/asm/clean/G3825.1a.fa -c G3825.1a -B '0 18834' \
  #                    -s realign $root/aln/raw/G3825.1a.bam $root/aln/filt/G3825.1a.bam \
  #                    $root/aln/tmp/G3825.1a
  paths['margin'] = 600
  if args.begin <= 2 and args.end >= 2:
    runner.run('bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {refname} -m {margin} '
               '-s realign -r {ref} {bam1raw} {bam1filt}'.format(**paths))
  else:
    print 'Skipping step 2.'

  # Remove duplicates
  #   samtools
  if args.begin <= 3 and args.end >= 3:
    dedup(paths['bam1filt'], paths['bam1dedup'], runner)
  else:
    print 'Skipping step 3.'

  # Extract reads
  #   Picard SamToFastq
  if args.begin <= 4 and args.end >= 4:
    runner.run('java -jar {picardir}/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT={bam1dedup} '
               'FASTQ={cleanfq1} SECOND_END_FASTQ={cleanfq2}'.format(**paths), ignore_err=True)
  else:
    print 'Skipping step 4.'

  # Assemble
  #   SPAdes
  #TODO: Set k-mers based on read length.
  #TODO: Check if successful (produces a contigs.fasta)
  # Example: $ spades.py --careful -k 21,33,55,77 -1 $FASTQ_DIR/${sample}_1.fastq \
  #            -2 $FASTQ_DIR/${sample}_2.fastq -o $ROOT/asm/orig/$sample
  if args.begin <= 5 and args.end >= 5:
    runner.run('spades.py --careful -k 21,33,55,77 -1 {cleanfq1} -2 {cleanfq2} -o {asmdir}'
               .format(**paths))
  else:
    print 'Skipping step 5.'

  # Clean assembly
  #   asm-unifier.py
  if args.begin <= 6 and args.end >= 6:
    runner.run('asm-unifier.py -n {sample} {ref} {asmdir}/contigs.fasta -o {asm} -l {asmlog}'
               .format(**paths))
  else:
    print 'Skipping step 6.'

  # Align assembly to reference
  #   LASTZ
  if args.begin <= 7 and args.end >= 7:
    runner.run('lastz {ref} {asm} > {lav}'.format(**paths))
  else:
    print 'Skipping step 7.'

  # Align to assembly
  #   BWA-MEM
  if args.begin <= 8 and args.end >= 8:
    align(paths['cleanfq1'], paths['cleanfq2'], paths['asm'], paths['bam2raw'], args.sample, runner)
  else:
    print 'Skipping step 8.'

  # Filter alignment
  #   pre-process-mt.sh
  # See notes for step 2.
  # In this second time around, might want to omit -s realign (do NM-edits filtering too).
  if args.begin <= 9 and args.end >= 9:
    runner.run('bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {sample} -m {margin} '
               '-s realign -r {asm} {bam2raw} {bam2filt}'.format(**paths))
  else:
    print 'Skipping step 9.'

  # Remove duplicates
  #   samtools
  if args.begin <= 10 and args.end >= 10:
    dedup(paths['bam2filt'], paths['bam2dedup'], runner)
  else:
    print 'Skipping step 10.'

  # Naive Variant Caller
  #TODO: Allow setting -q and -m.
  if args.begin <= 11 and args.end >= 11:
    runner.run('naive_variant_caller.py -q 30 -m 20 --ploidy 2 --use_strand '
               '--coverage_dtype uint32 --allow_out_of_bounds_positions --bam {bam2dedup} '
               '--index {bam2dedup}.bai -r {asm} -o {nvc}'.format(**paths))
  else:
    print 'Skipping step 11.'

  # nvc-filter.py
  #TODO: Allow changing coverage and frequency threshold.
  #TODO: Allow for samples with no indels (empty {nvc}).
  paths['cvg'] = 1000
  paths['freq'] = 1
  if args.begin <= 12 and args.end >= 12:
    runner.run('nvc-filter.py -r S -c {cvg} -f {freq} {nvc} > {nvcfilt}'.format(**paths))
  else:
    print 'Skipping step 12.'

  # inspect-reads.py
  # Note: -S overwrites read group names with the given sample name. Don't use in multi-sample
  # analysis.
  #TODO: Allow changing strand and mate bias thresholds.
  paths['strand'] = 1
  paths['mate'] = 1
  if args.begin <= 13 and args.end >= 13:
    runner.run('inspect-reads.py -tl -s {strand} -m {mate} -S {sample} {bam2dedup} -V {nvcfilt} '
               '-r {asm} > {vars1asm}'.format(**paths))
  else:
    print 'Skipping step 13.'

  # quick-liftover.py
  if args.begin <= 14 and args.end >= 14:
    runner.run('quick-liftover.py {lav} {vars1asm} > {vars1}'.format(**paths))
  else:
    print 'Skipping step 14.'

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
  paths = {'rg':rg_line, 'ref':ref, 'fq1':fastq1, 'fq2':fastq2, 'outbam':outbam, 'base':base}
  # Index reference, if needed.
  for ext in ('.amb', '.ann', '.bwt', '.sa', '.pac'):
    if not os.path.isfile(ref+ext):
      #TODO: determine algorithm based on file size ("is" if size is less than 2000000000 bytes).
      algorithm = 'bwtsw'
      runner.run('bwa index -a {algo} {ref}'.format(algo=algorithm, ref=ref))
      break
  runner.run('bwa mem -M -t 16 -R {rg} {ref} {fq1} {fq2} > {base}.sam'.format(**paths))
  runner.run('samtools view -Sb {base}.sam > {base}.tmp.bam'.format(**paths))
  runner.run('samtools sort {base}.tmp.bam {base}'.format(**paths))
  runner.run('samtools index {base}.bam'.format(**paths))
  runner.run('rm {base}.sam {base}.tmp.bam'.format(**paths))
  return outbam


def dedup(inbam, outbam, runner):
  (base, ext) = os.path.splitext(outbam)
  paths = {'inbam':inbam, 'base':base}
  runner.run('samtools view -b -F 1024 {inbam} > {base}.tmp.bam'.format(**paths))
  runner.run('samtools sort {base}.tmp.bam {base}'.format(**paths))
  runner.run('samtools index {base}.bam'.format(**paths))
  runner.run('rm {base}.tmp.bam'.format(**paths))


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
      except subprocess.CalledProcessError:
        if not ignore_err:
          raise
        else:
          sys.stderr.write('non-zero exit status on command\n$ '+command+'\n')
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
