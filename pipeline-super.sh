#!/usr/bin/env bash
set -ue

CHROM_DEFAULT="chrM"
BOUNDS_DEFAULT="600 16000"
# must be in PATH
REQUIRED_COMMANDS="awk bwa samtools"
# must be in same directory as this script
REQUIRED_SCRIPTS="pre-process-mt.sh"
PLATFORM=${PLATFORM:="ILLUMINA"}
DEBUG=${DEBUG:=}

USAGE="Usage: \$ $(basename $0) [options] ref.fa reads_1.fq reads_2.fq
Options:
-s sample:  Give the sample name, instead of inferring it from the first fastq
            filename.
-d dirname: The output directory to put the results (and intermediate files).
            The directory must already exist, and be empty.
-c chrom:   A chromosome to target the analysis to. Default: $CHROM_DEFAULT.
-B \"upper lower\": Bounds to hand to rm_chim_in_pair.py. Necessary, if using -c.
            Default: $BOUNDS_DEFAULT"

#TODO: find actual location of script, resolving links
scriptdir=$(dirname $0)

function main {

  # Check for required commands and scripts
  for cmd in $REQUIRED_COMMANDS; do
    if ! which $cmd >/dev/null 2>/dev/null; then
      fail "Error: could not find $cmd in PATH"
    fi
  done
  for script in $REQUIRED_SCRIPTS; do
    if [[ ! -f "$scriptdir/$script" ]]; then
      fail "Error: Required script $script missing from local directory."
    fi
  done

  # get options
  dir=''
  sample=''
  chrom="$CHROM_DEFAULT"
  chim_bounds="$BOUNDS_DEFAULT"
  while getopts ":s:d:c:B:h" opt; do
    case "$opt" in
      s) sample="$OPTARG";;
      d) dir="$OPTARG";;
      c) chrom="$OPTARG";;
      B) chim_bounds="$OPTARG";;
      h) fail "$USAGE";;
    esac
  done

  # get positional arguments
  narg=$OPTIND
  while [[ $narg -le $# ]]; do
    arg=${@:$narg:1}
    if [[ ${arg:0:1} == '-' ]]; then
      fail "Error: options like $arg must come before positional arguments."
    fi
    narg=$((narg+1))
  done
  positionals=$((narg-OPTIND))
  if [[ $positionals -lt 3 ]]; then
    fail "$USAGE"
  fi
  ref="${@:$OPTIND:1}"
  fastq1="${@:$OPTIND+1:1}"
  fastq2="${@:$OPTIND+2:1}"

  # If no -s, try to derive a sample name from the first fastq file.
  if [[ ! $sample ]]; then
    sample=$fastq1
    sample=$(basename $fastq1 .gz)
    sample=$(basename $sample .fq)
    sample=$(basename $sample .fastq)
    sample=$(basename $sample _1)
    sample=$(basename $sample _2)
    sample=$(echo "$sample" | sed -E 's/_S[0-9]+_L001_R[12]_001//')
  fi

  # Check that input files exist
  #TODO: Check that reference is indexed (avoid unclear BWA error message).
  for file in $ref $fastq1 $fastq2; do
    if [[ ! -s $file ]]; then
      fail "Error: \"$file\" nonexistent, inaccessible, or empty."
    fi
  done

  if [[ ! $dir ]]; then
    dir=$(dirname $fastq1)/$sample
  fi
  if [[ -e $dir ]]; then
    if [[ ! -d $dir ]] || [[ $(ls $dir) ]]; then
      fail "Error: \"$dir\" is either not a directory or not empty."
    fi
  else
    exho "mkdir $dir"
  fi

  # determine filenames
  status="$dir/status.txt"
  raw="$dir/raw.bam"
  filt="$dir/filt.bam"

  echo "\
scripts: $scriptdir
outdir:  $dir
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
sample:  $sample"

  echo -e "start\t$(date +%s)\t$(date)" > $status

  # Align to reference
  map $fastq1 $fastq2 $ref $raw $sample

  # Filter alignment
  exho "bash $scriptdir/pre-process-mt.sh -r $ref -c $chrom -B \"$chim_bounds\" $raw $filt"

  echo -e "end\t$(date +%s)\t$(date)" > $status

}


function map {
  fastq1=$1
  fastq2=$2
  ref=$3
  bam=$4
  sample=$5
  dir=$(dirname $bam)
  base=$(basename $bam .bam)
  exho "bwa mem -M -t 16 -R '@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM' $ref $fastq1 $fastq2 > $dir/$base.sam"
  exho "samtools view -Sb $dir/$base.sam > $dir/$base.tmp.bam"
  exho "samtools sort $dir/$base.tmp.bam $dir/$base"
  exho "samtools index $bam"
  exho "rm $dir/$base.sam $dir/$base.tmp.bam"
}

function fail {
  echo "$@" >&2
  exit 1
}

function exho {
  echo "+ \$ $1"
  if [[ ! $DEBUG ]]; then
    eval "$1"
  fi
}

main "$@"
