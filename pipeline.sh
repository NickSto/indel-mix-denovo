#!/usr/bin/env bash
if [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

CHROM_DEFAULT="chrM"
# must be in PATH
REQUIRED_COMMANDS=""
# must be in same directory as this script
REQUIRED_SCRIPTS=""
DEBUG=${DEBUG:=}

USAGE="Usage: \$ $(basename $0) alignment.bam"

#TODO: find actual location of script, resolving links
scriptdir=$(dirname $0)

function main {

  #################### SETUP ####################

  # Check for required commands and scripts.
  check_required "$scriptdir"

  # Parse command line arguments.
  read aln chrom <<< $(getmyopts "$@")

  # Try to determine sample name, if none given.
  # sample=$(get_sample "$sample" "$fastq1")

  # Make sure that input files and output directories exist.
  # dir=$(check_paths "$ref" "$fastq1" "$fastq2" "$dir")

  echo "\
scripts: $scriptdir
outdir:  $dir
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
sample:  $sample"


  #################### PIPELINE ####################

  echo -e "start\t$(date +%s)\t$(date)" > $status

  # Process BAM with NVC
  # exho "naive_variant_caller.py -q 30 -m 20 --ploidy 2 --use_strand \
  #   --coverage_dtype uint32 --allow_out_of_bounds_positions -r $ref \
  #   --region $chrom --bam $asm_align_filt --index $asm_align_filt.bai -o $vcf"

  # Process NVC output, filter for indels above $thres
  # exho "nvc-filter.py -r S -c 1000 -f 0.75 $nvcout > $root/indels/nvc/$family-filt.vcf"

  # Get read statistics
  # exho "inspect-reads.py -tl -S $sample $sampdir/$sample.bam -V $root/indels/nvc/$family-filt.vcf -r $asm > $root/indels/vars/$family/$sample-unfilt-asm.tsv"

  # Convert coordinates
  # exho "quick-liftover.py $lav $root/indels/vars/$family/$sample-unfilt-asm.tsv > $root/indels/vars/$family/$sample-unfilt.tsv"

  # Filter for strand and mate bias
  # exho "group-filter.py -s 1 -m 1 $sample_vars"

  echo -e "end\t$(date +%s)\t$(date)" >> $status

}


# Read in command line arguments
function getmyopts {
  # get options
  ## Default values can't be empty, otherwise it will be skipped when reading
  ## output of this function.
  dir='__NONE__'
  sample='__NONE__'
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
  aln="${@:$OPTIND:1}"
  # fastq1="${@:$OPTIND+1:1}"
  # fastq2="${@:$OPTIND+2:1}"
  echo "$aln" "$chrom"
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
    sample=$(echo "$sample" | sed -E 's/_S[0-9]+_L001_R[12]_001//')
  fi
  echo "$sample"
}


# Make sure that input files and output directories exist.
function check_paths {
  read ref fastq1 fastq2 dir <<< "$@"
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
  echo "$dir"
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
