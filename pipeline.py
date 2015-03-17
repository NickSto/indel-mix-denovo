#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess
import distutils.spawn

PLATFORM = 'ILLUMINA'
PICARDIR_DEFAULT = 'src/picard-tools-1.100'
REQUIRED_COMMANDS = ('rm', 'bwa', 'samtools', 'spades.py', 'lastz', 'naive_variant_caller.py')
REQUIRED_SCRIPTS = ('heteroplasmy/pre-process-mt.sh', 'asm-unifier.py', 'nvc-filter.py',
                    'inspect-reads.py', 'quick-liftover.py')
REQUIRED_PICARDS = ('SamToFastq.jar',)

OPT_DEFAULTS = {'begin':0, 'end':14, 'freq':1, 'cvg':1000, 'strand':1, 'mate':1,
                'margin':600, 'kmers':'21,33,55,77', 'read_length_minimum':40}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
FILENAMES = ('bam1raw.bam', 'bam1filt.bam', 'bam1dedup.bam', 'cleanfq1.fq', 'cleanfq2.fq',
             'asmdir', 'asm.fa', 'asmlog.log', 'lav.lav', 'bam2raw.bam', 'bam2filt.bam',
             'bam2dedup.bam', 'nvc.vcf', 'nvcfilt.vcf', 'vars_asm.tsv', 'vars.tsv')


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
    help='The sample id.')
  parser.add_argument('-r', '--refname', metavar='chrName', required=True,
    help='Id of reference sequence. Reads will be filtered for those that map to this sequence.')
  parser.add_argument('-R', '--filter-ref', metavar='filter-ref.fa',
    help='A reference consisting of the reference genome plus any other mapping targets you want '
         'to include for filtering purposes. Will filter the reads in the first step by mapping '
         'to this reference, then selecting only reads which map to the sequence identified by '
         '--refname. Optional. Will use the main reference genome by default.')
  parser.add_argument('-N', '--simulate', action='store_true',
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
    help='Read length. Default: "%(default)s"')
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
  params['scriptdir'] = os.path.relpath(os.path.dirname(os.path.realpath(sys.argv[0])))
  params['picardir'] = os.path.join(os.path.expanduser('~'), PICARDIR_DEFAULT)
  # Add arguments from command line.
  for arg in dir(args):
    if arg.startswith('_'):
      continue
    assert arg not in params, '{} in params. value: {}'.format(arg, params[arg])
    params[arg] = getattr(args, arg)
  # Compute some special arguments.
  if not args.filter_ref:
    params['filter_ref'] = args.ref
  if not params.get('rlen'):
    #TODO: Detect read length of FASTA files.
    raise Exception('Need to provide --read-length.')
  params['rlen_thres'] = int(round(params['rlen'] * params['read_length_minimum'] / 100))

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


  #################### PIPELINE STEPS ####################

  #TODO: Maybe eliminate some repetition and store command format strings in a simple list?

  step = 'Step 1: Align to reference with BWA MEM'
  if args.begin <= 1 and args.end >= 1:
    print step+':'
    align(args.fastq1, args.fastq2, params['filter_ref'], params['bam1raw'], args.sample, runner)
  elif args.end > 1:
    print 'Skipping '+step+'.'

  step = 'Step 2: Filter alignment with pre-process-mt.sh'
  #TODO: only use -s realign when necessary
  # Example command: $ pre-process-mt.sh -r $root/asm/clean/G3825.1a.fa -c G3825.1a -B '0 18834' \
  #                    -s realign $root/aln/raw/G3825.1a.bam $root/aln/filt/G3825.1a.bam \
  #                    $root/aln/tmp/G3825.1a
  if args.begin <= 2 and args.end >= 2:
    print step+':'
    runner.run('bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {refname} -m {margin} '
               '-s realign -r {ref} {bam1raw} {bam1filt}'.format(**params))
  elif args.end > 2:
    print 'Skipping '+step+'.'

  step = 'Step 3: Remove duplicates with Samtools'
  if args.begin <= 3 and args.end >= 3:
    print step+':'
    dedup(params['bam1filt'], params['bam1dedup'], runner)
  elif args.end > 3:
    print 'Skipping '+step+'.'

  step = 'Step 4: Extract reads with Picard SamToFastq'
  if args.begin <= 4 and args.end >= 4:
    print step+':'
    runner.run('java -jar {picardir}/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT={bam1dedup} '
               'FASTQ={cleanfq1} SECOND_END_FASTQ={cleanfq2}'.format(**params), ignore_err=True)
  elif args.end > 4:
    print 'Skipping '+step+'.'

  step = 'Step 5: Assemble with SPAdes'
  #TODO: Set k-mers based on read length.
  #TODO: Check if successful (produces a contigs.fasta)
  # Example: $ spades.py --careful -k 21,33,55,77 -1 $FASTQ_DIR/${sample}_1.fastq \
  #            -2 $FASTQ_DIR/${sample}_2.fastq -o $ROOT/asm/orig/$sample
  if args.begin <= 5 and args.end >= 5:
    print step+':'
    runner.run('spades.py --careful -k {kmers} -1 {cleanfq1} -2 {cleanfq2} -o {asmdir}'
               .format(**params))
  elif args.end > 5:
    print 'Skipping '+step+'.'

  step = 'Step 6: Clean assembly with asm-unifier.py'
  if args.begin <= 6 and args.end >= 6:
    print step+':'
    runner.run('python {scriptdir}/asm-unifier.py -n {sample} {ref} {asmdir}/contigs.fasta '
               '-o {asm} -l {asmlog}'.format(**params))
  elif args.end > 6:
    print 'Skipping '+step+'.'

  step = 'Step 7: Align assembly to reference with LASTZ'
  if args.begin <= 7 and args.end >= 7:
    print step+':'
    runner.run('lastz {ref} {asm} > {lav}'.format(**params))
  elif args.end > 7:
    print 'Skipping '+step+'.'

  step = 'Step 8: Align to assembly with BWA MEM'
  if args.begin <= 8 and args.end >= 8:
    print step+':'
    align(params['cleanfq1'], params['cleanfq2'], params['asm'], params['bam2raw'], args.sample,
          runner)
  elif args.end > 8:
    print 'Skipping '+step+'.'

  step = 'Step 9: Filter alignment with pre-process-mt.sh'
  # See notes for step 2.
  # In this second time around, might want to omit -s realign (do NM-edits filtering too).
  if args.begin <= 9 and args.end >= 9:
    print step+':'
    runner.run('bash {scriptdir}/heteroplasmy/pre-process-mt.sh -c {sample} -m {margin} '
               '-s realign -r {asm} {bam2raw} {bam2filt}'.format(**params))
  elif args.end > 9:
    print 'Skipping '+step+'.'

  step = 'Step 10: Remove duplicates with Samtools:'
  if args.begin <= 10 and args.end >= 10:
    print step+':'
    dedup(params['bam2filt'], params['bam2dedup'], runner)
  elif args.end > 10:
    print 'Skipping '+step+'.'

  step = 'Step 11: Extract variants with Naive Variant Caller'
  #TODO: Allow setting -q and -m.
  if args.begin <= 11 and args.end >= 11:
    print step+':'
    runner.run('naive_variant_caller.py -q 30 -m 20 --ploidy 2 --use_strand '
               '--coverage_dtype uint32 --allow_out_of_bounds_positions --bam {bam2dedup} '
               '--index {bam2dedup}.bai -r {asm} -o {nvc}'.format(**params))
  elif args.end > 11:
    print 'Skipping '+step+'.'

  step = 'Step 12: Filter variants with nvc-filter.py'
  #TODO: Allow for samples with no indels (empty {nvc}).
  if args.begin <= 12 and args.end >= 12:
    print step+':'
    runner.run('python {scriptdir}/nvc-filter.py -r S -c {cvg} -f {freq} {nvc} > {nvcfilt}'
               .format(**params))
  elif args.end > 12:
    print 'Skipping '+step+'.'

  step = 'Step 13: Filter further by BAM stats with inspect-reads.py'
  # Note: -S overwrites read group names with the given sample name. Don't use in multi-sample
  # analysis.
  if args.begin <= 13 and args.end >= 13:
    print step+':'
    runner.run('python {scriptdir}/inspect-reads.py -tl -s {strand} -m {mate} -S {sample} '
               '{bam2dedup} -V {nvcfilt} -r {asm} > {vars_asm}'.format(**params))
  elif args.end > 13:
    print 'Skipping '+step+'.'

  step = 'Step 14: Convert to reference coordinates with quick-liftover.py'
  if args.begin <= 14 and args.end >= 14:
    print step+':'
    runner.run('python {scriptdir}/quick-liftover.py {lav} {vars_asm} > {vars}'.format(**params))
  elif args.end > 14:
    print 'Skipping '+step+'.'

  if args.end >= 14:
    print 'Finished. Final variants are in {vars}.'.format(**params)


def align(fastq1, fastq2, ref, outbam, sample, runner):
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
  runner.run('samtools index {base}.bam'.format(**params))
  runner.run('rm {base}.sam {base}.tmp.bam'.format(**params))
  return outbam


def dedup(inbam, outbam, runner):
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
