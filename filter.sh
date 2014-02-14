#!/usr/bin/env bash
set -ue

if [[ $# -lt 4 ]]; then
  echo "USAGE: $(basename $0) families.tsv choices.tsv family.bam nvc-out.vcf [asm.lav]
There must be a directory structure like:
root
 |- asm
 |   \- lastz
 \- indels
    |- filt - family.bam
    |- nvc  - nvc-out.vcf
    \- vars"
  exit 1
fi

families=$1
choices=$2
fambam=$3
nvcout=$4

function fail {
  echo $1 >&2
  exit 1
}


for command in samtools nvc-filter.py inspect-reads.py quick-liftover.py group-filter.py; do
  if ! which $command >/dev/null; then
    fail "Error: cannot find $command command."
  fi
done


##### check inputs #####

# family.bam must end in .bam
if [[ ${fambam:(-4)} != '.bam' ]]; then
  fail "Error: $fambam does not end in .bam"
fi
root=$(dirname $(dirname $(dirname $(readlink -e $fambam))))
# required directories
if [[ $(dirname $(readlink -e $fambam)) != $root/indels/filt ]]; then
  fail "Error in directory structure: dirname $fambam != $root/indels/filt"
fi
if [[ $(dirname $(readlink -e $nvcout)) != $root/indels/nvc ]]; then
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
family=$(basename $fambam .bam)
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
  samtools view -b -r $sample $fambam > $tmpdir/$sample.bam
done

# filter for indels above 0.75%
nvc-filter.py -r S -c 1000 -f 0.75 $nvcout > $root/indels/nvc/$family-filt.vcf

# get read statistics
for sample in $samples; do
  inspect-reads.py -t $tmpdir/$sample.bam -V $root/indels/nvc/$family-filt.vcf > $root/indels/vars/$family/$sample-unfilt-asm.tsv
done

# convert coordinates
sample_vars=''
for sample in $samples; do
  quick-liftover.py $lav $root/indels/vars/$family/$sample-asm-unfilt.tsv > $root/indels/vars/$family/$sample-unfilt.tsv
  sample_vars="$sample_vars $root/indels/vars/$family/$sample-unfilt.tsv"
done

# filter for strand bias and mate bias
group-filter.py -s 1 -m 1 $sample_vars

# rename output files
for sample in $samples; do
  mv $root/indels/vars/$family/$sample-unfilt-filt.tsv $root/indels/vars/$family/$sample.tsv
done
