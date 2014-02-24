#!/usr/bin/env bash
set -ue

# assumes form M477-filt1-filt2.bam or M477.bam
BAMREGEX='s/^([^-]+)(-.+)*\.bam$/\1/g'

function fail {
  echo $1 >&2
  exit 1
}

if [[ $# -lt 4 ]]; then
  fail "USAGE: $(basename $0) families.tsv choices.tsv family.bam nvc-out.vcf [asm.lav]
There must be a directory structure like:
root
 |- asm
 |   \- lastz
 \- indels
    |- filt - family.bam
    |- nvc  - nvc-out.vcf
    \- vars"
fi

families=$1
choices=$2
fambam=$3
nvcout=$4


for command in samtools nvc-filter.py inspect-reads.py quick-liftover.py group-filter.py; do
  if ! which $command >/dev/null; then
    fail "Error: cannot find $command command."
  fi
done


##### check inputs #####

root=$(dirname $(dirname $(cd $(dirname $fambam); pwd)))
# required directories
if [[ $(cd $(dirname $fambam); pwd) != $root/indels/filt ]]; then
  fail "Error in directory structure: dirname $fambam != $root/indels/filt"
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
# get asm file
if [[ $# -ge 5 ]]; then
  lav=$5
else
  asmname=$(grep -E "^$family\b" $choices | cut -f 2)
  if [[ ! "$asmname" ]]; then
    fail "Error: $family not found in $choices"
  fi
  lav=$root/asm/lastz/$(basename $asmname .fa)'.lav'
fi
if [[ ! -f $lav ]]; then
  fail "Error: cannot find file $lav"
fi
tmpdir=$root/indels/filt/tmp


##### do actual analysis #####

if [[ ! -d $tmpdir ]]; then
  mkdir $tmpdir
fi
if [[ ! -d $root/indels/vars/$family ]]; then
  mkdir -p $root/indels/vars/$family
fi

set -x

# break family BAM into individual sample BAMs
for sample in $samples; do
  if [[ ! -s $tmpdir/$sample.bam ]]; then
    samtools view -b -r $sample $fambam > $tmpdir/$sample.bam
  fi
done

# filter for indels above 0.75%
if [[ ! -s $root/indels/nvc/$family-filt.vcf ]]; then
  nvc-filter.py -r S -c 1000 -f 0.75 $nvcout > $root/indels/nvc/$family-filt.vcf
fi

# get read statistics
for sample in $samples; do
  if [[ ! -s $root/indels/vars/$family/$sample-unfilt-asm.tsv ]]; then
    inspect-reads.py -tl -S $sample $tmpdir/$sample.bam -V $root/indels/nvc/$family-filt.vcf > $root/indels/vars/$family/$sample-unfilt-asm.tsv
  fi
done

# convert coordinates
sample_vars=''
for sample in $samples; do
  if [[ ! -s $root/indels/vars/$family/$sample-unfilt.tsv ]]; then
    quick-liftover.py $lav $root/indels/vars/$family/$sample-unfilt-asm.tsv > $root/indels/vars/$family/$sample-unfilt.tsv
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
  group-filter.py -s 1 -m 1 $sample_vars
fi

# rename output files
for sample in $samples; do
  mv $root/indels/vars/$family/$sample-unfilt-filt.tsv $root/indels/vars/$family/$sample.tsv
done
