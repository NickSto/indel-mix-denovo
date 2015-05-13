#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

SAMPLE_DEFAULT='G3825.1a'

ROOT=/nfs/brubeck.bx.psu.edu/scratch4/nick/ebola
FASTQ_DIR=/galaxy/home/boris/boris4/ebola/fastqFiles
REF=/nfs/brubeck.bx.psu.edu/scratch4/nick/refs/simple/ebola/EBOV.fa

USAGE="Usage: \$ $(basename $0) [sampleid]
Execute simple series of steps to run the indel pipeline on a dataset. Don't run this without
reading the commands it'll execute! If an argument is given, it is used as the \$sample variable
instead of \$SAMPLE_DEFAULT.
Current settings:
\$SAMPLE_DEFAULT: $SAMPLE_DEFAULT
\$ROOT:           $ROOT
\$FASTQ_DIR:      $FASTQ_DIR
\$REF:            $REF"

function main {

  if [[ $# -gt 0 ]]; then
    if [[ $1 == '-h' ]]; then
      fail "$USAGE"
    fi
    sample=$1
  else
    sample=$SAMPLE_DEFAULT
  fi

  # map to reference
  map $FASTQ_DIR/${sample}_{1,2}.fastq $REF $ROOT/map/$sample/raw.bam $sample

  # assemble
  if ! [[ -f $FASTQ_DIR/${sample}_1.fastq ]] || ! [[ $FASTQ_DIR/${sample}_2.fastq ]]; then
    fail "Error: input fastq files not found."
  fi
  if [[ -e $ROOT/asm/orig/$sample ]] || ! [[ -d $ROOT/asm/orig ]]; then
    fail "Error: output assembly directory exists or its parent does not."
  fi
  set -x
  srun -C new spades.py --careful -k 21,33,55,77 -1 $FASTQ_DIR/${sample}_1.fastq \
    -2 $FASTQ_DIR/${sample}_2.fastq -o $ROOT/asm/orig/$sample &
  set - -x

  # clean assembly
  if ! [[ -f $REF ]] || ! [[ -f orig/$sample/contigs.fasta ]]; then
    fail "Error: assembly unifier input files not found."
  fi
  if [[ -e $ROOT/asm/clean/$sample.fa ]] || ! [[ -d $ROOT/asm/clean ]]; then
    fail "Error: assembly unifier output file already exists or its parent does not."
  fi
  set -x
  srun -C new asm-unifier.py -n $sample $REF $ROOT/asm/orig/$sample/contigs.fasta > \
    $ROOT/asm/clean/$sample.fa
  set - -x



}


# Map fastq's to reference with BWA
function map {
  read fastq1 fastq2 ref bam sample <<< "$@"
  dir=$(dirname $bam)
  base=$(basename $bam .bam)
  if [[ ! -s $ref.amb ]] || [[ ! -s $ref.ann ]] || [[ ! -s $ref.bwt ]] || \
      [[ ! -s $ref.sa ]] || [[ ! -s $ref.pac ]]; then
    algo="bwtsw"
    if [[ $(du -sb $ref | awk '{print $1}') -lt 2000000000 ]]; then
      algo="is"
    fi
    exho "bwa index -a $algo $ref"
  fi
  set -x
  bwa mem -M -t 16 -R '@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM' $ref $fastq1 $fastq2 > $dir/$base.sam
  samtools view -Sb $dir/$base.sam > $dir/$base.tmp.bam
  samtools sort $dir/$base.tmp.bam $dir/$base
  samtools index $bam
  rm $dir/$base.sam $dir/$base.tmp.bam
  set - -x
}


function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
