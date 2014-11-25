#!/usr/bin/env bash
if [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

CHROM_DEFAULT="chrM"
BOUNDS_DEFAULT="600 16000"
# must be in PATH
REQUIRED_COMMANDS="awk java samtools bwa lastz spades.py"
# must be in same directory as this script
REQUIRED_SCRIPTS="pre-process-mt.sh asm-unifier.py"
PICARD_DIR=${PICARD_DIR:-~/src/picard-tools-1.100}
PLATFORM=${PLATFORM:="ILLUMINA"}
DEBUG=${DEBUG:=}

USAGE="Usage: \$ $(basename $0) [options] ref.fa filter-ref.fa reads_1.fq reads_2.fq
Run the tasks required to produce the datasets the indel pipeline can run on.
This will filter, assemble, map, and filter the data from one sample. Then,
pipeline.sh will run the indel pipeline on the assembly-aligned data.
pipeline-meta.sh can be used to run this on a set of related samples, merge them
into a single BAM marked by read groups, then run the indel pipeline on it.
  Positional Arguments:
ref.fa:           The target reference sequence. This will be used to guide
                  cleanup of the de novo assembly.
filter-ref.fa:    A reference, composed of the target sequence, plus all other
                  sequences present in the sample which could potentially be a
                  source of reads. This will be the mapping target for reads in
                  the filtering step.
reads_1.fq:       Input reads, mate 1, in fastq format.
reads_2.fq:       Ditto, mate 2.
  Options:
-s sample:        Give the sample name, instead of inferring it from the first
                  fastq filename.
-d dirname:       The output directory to put the results (and intermediate
                  files). If the directory already exists, it must be empty.
-o out.bam:       Where to put/name the final, assembly-aligned BAM.
-c chrom:         A chromosome id to target the analysis to. Reads will be
                  filtered to those which align to this chromosome instead of
                  the rest of the filter-ref.fa. Default: $CHROM_DEFAULT.
-B \"upper lower\": Bounds to hand to rm_chim_in_pair.py. Necessary, if using -c.
                  Default: $BOUNDS_DEFAULT
N.B.: None of the arguments (except -B) can contain spaces."

#TODO: find actual location of script, resolving links
scriptdir=$(dirname $0)

function main {

  #################### SETUP ####################

  # Check for required commands and scripts.
  check_required "$scriptdir"

  # Parse command line arguments.
  read ref filter_ref fastq1 fastq2 dir sample outpath chrom bound1 bound2 \
    <<< $(getmyopts "$@")
  if [[ ! $ref ]]; then
    exit 1
  fi

  # Try to determine sample name, if none given.
  sample=$(get_sample $sample $fastq1)

  # Make sure that input files and output directories exist.
  dir=$(check_paths $ref $fastq1 $fastq2 $dir $sample)

  echo "\
scripts: $scriptdir
outdir:  $dir
outpath: $outpath
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
sample:  $sample"


  #################### PIPELINE ####################

  echo -e "start\t$(date +%s)\t$(date)" > $dir/status-pre.txt

  # Align to reference
  map $fastq1 $fastq2 $filter_ref $dir/ref_raw.bam $sample

  # Filter alignment
  exho "bash $scriptdir/pre-process-mt.sh -r $filter_ref -s chimrlen -c $chrom \
    -B \"$bound1 $bound2\" $dir/ref_raw.bam $dir/ref_filt.bam"

  # Extract reads
  set +e  # Picard exits with spurious errors
  exho "java -jar $PICARD_DIR/SamToFastq.jar INPUT=$dir/ref_filt.bam \
    FASTQ=$dir/filt_1.fq SECOND_END_FASTQ=$dir/filt_2.fq \
    VALIDATION_STRINGENCY=SILENT"
  set -e

  # Assemble
  exho "spades.py -k 21,33,55,77,99,127 --careful -1 $dir/filt_1.fq \
    -2 $dir/filt_2.fq -o $dir/asm"

  # Clean up assembly
  exho "python $scriptdir/asm-unifier.py -n $chrom $ref \
    $dir/asm/contigs.fasta -o $dir/asm.fa"

  # Map to assembly
  map $dir/filt_1.fq $dir/filt_2.fq $dir/asm.fa $dir/asm_raw.bam $sample

  # Filter assembly alignment
  exho "bash $scriptdir/pre-process-mt.sh -r $dir/asm.fa -s realign -c $chrom \
    -B \"$bound1 $bound2\" $dir/asm_raw.bam $dir/asm_filt.bam"

  if [[ $outpath != '__NONE__' ]]; then
    mv $dir/asm_filt.bam $outpath
  fi

  echo -e "end\t$(date +%s)\t$(date)" >> $dir/status-pre.txt

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
  exho "bwa mem -M -t 16 -R '@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM' $ref $fastq1 $fastq2 > $dir/$base.sam"
  exho "samtools view -Sb $dir/$base.sam > $dir/$base.tmp.bam"
  exho "samtools sort $dir/$base.tmp.bam $dir/$base"
  exho "samtools index $bam"
  exho "rm $dir/$base.sam $dir/$base.tmp.bam"
}


# Read in command line arguments
function getmyopts {
  # get options
  ## Default values can't be empty, otherwise it will be skipped when reading
  ## output of this function.
  dir='__NONE__'
  sample='__NONE__'
  outpath='__NONE__'
  chrom="$CHROM_DEFAULT"
  chim_bounds="$BOUNDS_DEFAULT"
  while getopts ":s:d:o:c:B:h" opt; do
    case "$opt" in
      s) sample="$OPTARG";;
      d) dir="$OPTARG";;
      o) outpath="$OPTARG";;
      c) chrom="$OPTARG";;
      B) chim_bounds="$OPTARG";;
      h) fail "$USAGE";;
    esac
  done
  read bound1 bound2 <<< $(echo $chim_bounds)
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
  filter_ref="${@:$OPTIND+1:1}"
  fastq1="${@:$OPTIND+2:1}"
  fastq2="${@:$OPTIND+3:1}"
  # Catch arguments with spaces
  if spaces "$ref" || spaces "$filter_ref" || spaces "$fastq1" \
      || spaces "$fastq2" || spaces "$dir" || spaces "$sample" \
      || spaces "$outpath" || spaces "$chrom"; then
    fail "Error: Arguments must not contain spaces"
  fi
  echo $ref $filter_ref $fastq1 $fastq2 $dir $sample $outpath $chrom $bound1 $bound2
}


# Does the first argument contain whitespace?
function spaces {
  value="$1"
  if [[ "$value" == "$(echo "$value" | tr -d ' \n\r\t')" ]]; then
    false
  else
    true
  fi
}


# Check for required commands and scripts
function check_required {
  scriptdir="$1"
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
  if [[ ! -d $PICARD_DIR ]] || [[ -z $(ls $PICARD_DIR) ]]; then
    fail "Picard directory $PICARD_DIR missing or empty."
  fi
}


# If no -s, try to derive a sample name from the first fastq file.
function get_sample {
  read sample fastq1 <<< "$@"
  if [[ $sample == '__NONE__' ]]; then
    sample=$fastq1
    sample=$(basename $fastq1 .gz)
    sample=$(basename $sample .fq)
    sample=$(basename $sample .fastq)
    sample=$(basename $sample _1)
    sample=$(basename $sample _2)
    sample=$(echo $sample | sed -E 's/_S[0-9]+_L001_R[12]_001//')
  fi
  echo $sample
}


# Make sure that input files and output directories exist.
function check_paths {
  read ref fastq1 fastq2 dir sample <<< "$@"
  #TODO: Check that reference is indexed (avoid unclear BWA error message).
  # Check for existence of input files.
  for file in $ref $fastq1 $fastq2; do
    if [[ ! -s $file ]]; then
      fail "Error: \"$file\" nonexistent, inaccessible, or empty."
    fi
  done
  # Check $dir.
  # If not given with -d, determine it from the fastq1 path.
  # Make sure if it exists, it's empty.
  if [[ $dir == '__NONE__' ]]; then
    dir=$(dirname $fastq1)/$sample
  fi
  if [[ -e $dir ]]; then
    if [[ ! -d $dir ]] || [[ $(ls $dir) ]]; then
      fail "Error: \"$dir\" is either not a directory or not empty."
    fi
  else
    exho "mkdir -p $dir" >&2
  fi
  echo $dir
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
