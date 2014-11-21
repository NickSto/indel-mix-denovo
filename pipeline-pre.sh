#!/usr/bin/env bash
if [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

CHROM_DEFAULT="chrM"
BOUNDS_DEFAULT="600 16000"
# must be in PATH
REQUIRED_COMMANDS="awk bwa samtools java spades.py"
# must be in same directory as this script
REQUIRED_SCRIPTS="pre-process-mt.sh asm-unifier.py"
PICARD_DIR=${PICARD_DIR:-~/src/picard-tools-1.100}
PLATFORM=${PLATFORM:="ILLUMINA"}
DEBUG=${DEBUG:=}

USAGE="Usage: \$ $(basename $0) [options] ref.fa reads_1.fq reads_2.fq
Run the tasks required to produce the datasets the indel pipeline can run on.
This will filter, assemble, map, and filter the data from one sample. Then,
pipeline.sh will run the indel pipeline on the assembly-aligned data.
pipeline-meta.sh can be used to run this on a set of related samples, merge them
into a single BAM marked by read groups, then run the indel pipeline on it.
Options:
-s sample:  Give the sample name, instead of inferring it from the first fastq
            filename.
-d dirname: The output directory to put the results (and intermediate files).
            If the directory already exists, it must be empty.
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
  if [[ ! -d $PICARD_DIR ]] || [[ -z $(ls $PICARD_DIR) ]]; then
    fail "Picard directory $PICARD_DIR missing or empty."
  fi

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
    exho "mkdir -p $dir"
  fi

  echo "\
scripts: $scriptdir
outdir:  $dir
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
sample:  $sample"

  # determine filenames
  status="$dir/status.txt"
  ref_align_raw="$dir/ref_raw.bam"
  ref_align_filt="$dir/ref_filt.bam"
  fastq1_filt="$dir/filt_1.fq"
  fastq2_filt="$dir/filt_2.fq"
  asm_dir="$dir/asm"
  asm="$dir/asm.fa"
  asm_align_raw="$dir/asm_raw.bam"
  asm_align_filt="$dir/asm_filt.bam"

  echo -e "start\t$(date +%s)\t$(date)" > $status

  # Align to reference
  map $fastq1 $fastq2 $ref $ref_align_raw $sample

  # Filter alignment
  exho "bash $scriptdir/pre-process-mt.sh -r $ref -s chimrlen -c $chrom \
    -B \"$chim_bounds\" $ref_align_raw $ref_align_filt"

  # Extract reads
  set +e  # Picard exits with spurious errors
  exho "java -jar $PICARD_DIR/SamToFastq.jar INPUT=$ref_align_filt \
    FASTQ=$fastq1_filt SECOND_END_FASTQ=$fastq2_filt \
    VALIDATION_STRINGENCY=SILENT"
  set -e

  # Assemble
  exho "spades.py -k 21,33,55,77,99,127 --careful -1 $fastq1_filt \
    -2 $fastq2_filt -o $asm_dir"

  # Clean up assembly
  exho "python $scriptdir/asm-unifier.py -n $sample $ref \
    $asm_dir/contigs.fasta -o $asm"

  # Map to assembly
  map $fastq1_filt $fastq2_filt $asm $asm_align_raw $sample

  # Filter assembly alignment
  exho "bash $scriptdir/pre-process-mt.sh -r $ref -s realign -c $chrom -B \"$chim_bounds\" \
    $asm_align_raw $asm_align_filt"

  #TODO: After this, I need to merge multiple individuals into a family.bam,
  #      then run the actual indel pipeline on that.

  # Process BAM with NVC
  # exho "naive_variant_caller.py -q 30 -m 20 --ploidy 2 --use_strand \
  #   --coverage_dtype uint32 --allow_out_of_bounds_positions -r $ref \
  #   --region $chrom --bam $asm_align_filt --index $asm_align_filt.bai -o $vcf"

  # Process NVC output with Variant Annotator, filter for $chrom
  # exho "nvc-filter.py -r S -c 1000 -f 0.75 -n -i $vcf | awk '$AWK_FILT' > $vars"

  echo -e "end\t$(date +%s)\t$(date)" >> $status

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
