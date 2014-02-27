#!/usr/bin/env bash
set -e

# assumes form M477-filt1-filt2.bam or M477.bam
BAMREGEX='s/^([^-]+)(-.+)*\.bam$/\1/g'

function fail {
  echo $1 >&2
  exit 1
}
function ecec {
  echo "$1"
  exec "$1"
}

if [[ $# -lt 4 ]]; then
  echo "USAGE: $(basename $0) families.tsv choices.tsv family.bam nvc-out.vcf [asm.lav] [asm.fa]
There must be a directory structure like:
root
 |- asm
 |   |- curated - asm.fa
 |   \- lastz - asm.lav
 \- indels
    |- nvc  - nvc-out.vcf
    |- vars
    \- filt
        |- family - family.bam
        \- sample
The family BAMs will be decomposed into sample BAMs and stored in
root/indels/filt/sample. If a sample BAM of the right name already exists there,
it will use that instead." >&2
  exit 1
fi

families=$1
choices=$2
fambam=$3
nvcout=$4
lav=$5
asm=$6

for command in samtools nvc-filter.py inspect-reads.py quick-liftover.py group-filter.py; do
  if ! which $command >/dev/null; then
    fail "Error: cannot find $command command."
  fi
done


##### check inputs #####

root=$(dirname $(dirname $(dirname $(cd $(dirname $fambam); pwd))))
# required directories
if [[ $(cd $(dirname $fambam); pwd) != $root/indels/filt/family ]]; then
  fail "Error in directory structure: dirname $fambam != $root/indels/filt/family"
fi
if [[ $(cd $(dirname $nvcout); pwd) != $root/indels/nvc ]]; then
  fail "Error in directory structure: dirname $nvcout != $root/indels/nvc"
fi
# required files
for file in $families $choices $fambam $nvcout; do
  if [[ ! -f $file ]]; then
    fail "Error: cannot find file $file"
  fi
done


##### establish derived parameters #####

# get family and sample names
family=$(basename $fambam | sed -E $BAMREGEX)
if [[ ! "$family" ]]; then
  fail "Error: $fambam does not match pattern $BAMREGEX"
fi
samples=$(grep -E "^$family\b" $families | cut -f 2-)
if [[ ! "$samples" ]]; then
  fail "Error: $family not found in $families"
fi
echo "Working on family $family, samples $samples"
# get the base name of the assembly file
asmname=$(basename $(grep -E "^$family\b" $choices | cut -f 2) .fa)
if [[ -z $asmname ]] && ( [[ -z $lav ]] || [[ -z $asm ]] ); then
  fail "Error: $family not found in $choices"
fi
# get lav file if not provided already
if [[ -z $lav ]]; then
  lav=$root/asm/lastz/$asmname.lav
fi
if [[ ! -f $lav ]]; then
  fail "Error: cannot find file $lav"
fi
# get asm file if not provided already
if [[ -z $asm ]]; then
  asm=$root/asm/curated/$asmname.fa
fi
if [[ ! -f $asm ]]; then
  fail "Error: cannot find file $asm"
fi
sampdir=$root/indels/filt/sample


##### do actual analysis #####

if [[ ! -d $sampdir ]]; then
  ecec "mkdir $sampdir"
fi
if [[ ! -d $root/indels/vars/$family ]]; then
  ecec "mkdir -p $root/indels/vars/$family"
fi

# break family BAM into individual sample BAMs
for sample in $samples; do
  if [[ ! -s $sampdir/$sample.bam ]]; then
    ecec "samtools view -b -r $sample $fambam > $sampdir/$sample.bam"
  fi
done

# filter for indels above 0.75%
if [[ ! -s $root/indels/nvc/$family-filt.vcf ]]; then
  ecec "nvc-filter.py -r S -c 1000 -f 0.75 $nvcout > $root/indels/nvc/$family-filt.vcf"
fi

# get read statistics
for sample in $samples; do
  if [[ ! -s $root/indels/vars/$family/$sample-unfilt-asm.tsv ]]; then
    ecec "inspect-reads.py -tl -S $sample $sampdir/$sample.bam -V $root/indels/nvc/$family-filt.vcf -r $asm > $root/indels/vars/$family/$sample-unfilt-asm.tsv"
  fi
done

# convert coordinates
sample_vars=''
for sample in $samples; do
  if [[ ! -s $root/indels/vars/$family/$sample-unfilt.tsv ]]; then
    ecec "quick-liftover.py $lav $root/indels/vars/$family/$sample-unfilt-asm.tsv > $root/indels/vars/$family/$sample-unfilt.tsv"
  fi
  sample_vars="$sample_vars $root/indels/vars/$family/$sample-unfilt.tsv"
done

# filter for strand bias and mate bias
do_filter=""
for sample in $samples; do
  if [[ ! -s $root/indels/vars/$family/$sample.tsv ]]; then
    do_filter="true"
  fi
done
if [[ "$do_filter" ]]; then
  ecec "group-filter.py -s 1 -m 1 $sample_vars"
fi

# rename output files
for sample in $samples; do
  ecec "mv $root/indels/vars/$family/$sample-unfilt-filt.tsv $root/indels/vars/$family/$sample.tsv"
done
