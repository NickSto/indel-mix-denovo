#!/usr/bin/env bash
set -ue

# note: the reference does not have to be indexed; bamleftalign will do that.
REF_DEFAULT="$HOME/bx/data/chrM-rCRS.fa"
CHROM_DEFAULT="chrM"
BOUNDS_DEFAULT="600 16000"
REQUIRED_SCRIPTS="rm_chim_in_pair.py nm-ratio.select.py get_major_from_bam.py"
REQUIRED_COMMANDS="java bamtools samtools bamleftalign"
PICARD_DIR=${PICARD_DIR:-~/src/picard-tools-1.100}
REGEX="[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*."
USAGE="Usage: \$ $(basename $0) [-r reference.fa] input.bam [output-name.bam [tmp_dir_name]]
N.B.: If you give a temporary directory name, it will not be deleted at the end.
-r: The reference sequence. Default: $REF_DEFAULT
-s: A processing step to stop at. Will return that result file instead.
    Step names:
      filter1
      chimrlen
      realign
-c: A chromosome to target the analysis to. Default: $CHROM_DEFAULT.
    N.B.: Not fully implemented yet! Only works up through chimrlen!
-B: Bounds to hand to rm_chim_in_pair.py. Necessary, if using -c. Default:
    $BOUNDS_DEFAULT"

#TODO: find actual location of script, resolving links
scriptdir=$(dirname $0)

function fail {
  echo "Error: $1" >&2
  exit 1
}

# read in arguments
ref="$REF_DEFAULT"
chrom="$CHROM_DEFAULT"
chim_bounds="$BOUNDS_DEFAULT"
stopat=''
while getopts ":r:s:c:B:h" opt; do
  case "$opt" in
    r) ref="$OPTARG";;
    s) stopat="$OPTARG";;
    c) chrom="$OPTARG";;
    B) chim_bounds="$OPTARG";;
    h) echo "$USAGE" >&2
       exit 1;;
  esac
done
inpath="${@:$OPTIND:1}"     # input BAM: aligned to hg19 + rCRS + spike-ins
outfile="${@:$OPTIND+1:1}"
tmpdir="${@:$OPTIND+2:1}"
# validate arguments
  #TODO: check for valid $stopat
if [[ ! $inpath ]]; then
  echo "$USAGE" >&2
  exit 1
fi
bamdir=$(dirname "$inpath")
name=$(basename "$inpath" .bam)
if [[ ! $outfile ]]; then
  outfile="$bamdir/$name.filt.bam"
fi

# check for required commands, scripts, and files
for command in $REQUIRED_COMMANDS; do
  if ! which $command >/dev/null 2>/dev/null; then
    fail "Required command $command missing from PATH."
  fi
done
for script in $REQUIRED_SCRIPTS; do
  if [[ ! -f "$scriptdir/$script" ]]; then
    fail "Required script $script missing from local directory."
  fi
done
if [[ ! -d $PICARD_DIR ]] || [[ -z $(ls $PICARD_DIR) ]]; then
  fail "Picard directory $PICARD_DIR missing or empty."
fi
if [[ ! -s "$inpath" ]]; then
  fail "input BAM \"$inpath\" not found or empty"
fi
type=$(file -b --mime-type "$inpath")
if [[ $type != "application/gzip" ]] && [[ $type != "application/x-gzip" ]]; then
  fail "input file \"$inpath\" type is \"$type\", must be \"application/gzip\""
fi
if [[ ! -s "$ref" ]]; then
  fail "Reference file \"$ref\" missing."
fi

# make a tmp directory
if [[ $tmpdir ]]; then
  keep_tmp='true'
  if [[ -e $tmpdir ]]; then
    fail "temp directory \"$tmpdir\" already exists."
  fi
else
  keep_tmp=''
  tmpdir="$bamdir/$name-tmp"
  tries=0
  while [[ -e "$tmpdir" ]]; do
    tmpdir="$bamdir/$name-tmp-$RANDOM"
    if [[ $((tries++)) -gt 20 ]]; then
      fail "temp directory \"$tmpdir\" already exists."
    fi 
  done
fi
mkdir "$tmpdir"


# Save the given output file, remove the tmp dir, and exit.
# The first argument is the full path to the final output file (in the tmp dir).
# If it ends in .bam, any index file will also be moved and renamed.
function finish {
  tmp_outfile="$1"

  # move the final output file out of the temp folder and rename it
  mv "$tmp_outfile" "$outfile"
  # if it's a BAM, move its index too
  if [[ $(basename "$tmp_outfile" .bam) != "$tmp_outfile" ]] && [[ -f "$tmp_outfile.bai" ]]; then
    mv "$tmp_outfile.bai" "$outfile.bai"
  fi

  # cleanup
  echo "--- removing intermediate files ---"
  if [[ ! $keep_tmp ]]; then
    rm -r "$tmpdir"
  fi

  echo "--- done ---"
  exit
}


#################### PIPELINE ####################

echo "--- sort by coordinate started ---"
input="$inpath"
output="sorted1.bam"
java -jar $PICARD_DIR/SortSam.jar \
  I="$input"                      \
  O="$tmpdir/$output"             \
  SO=coordinate

# Mark duplicates
echo "--- markDup started ---"
input="sorted1.bam"
output="dupmarked.bam"
output2="dupmarked.metrics"
java -jar $PICARD_DIR/MarkDuplicates.jar \
  I="$tmpdir/$input"                     \
  O="$tmpdir/$output"                    \
  METRICS_FILE="$tmpdir/$output2"        \
  READ_NAME_REGEX="$REGEX"

# Filter for mapped in proper pair
echo "--- Selection of proper started ---"
input="dupmarked.bam"
output="filtered1.bam"
bamtools filter            \
  -region $chrom           \
  -isPaired true           \
  -isProperPair true       \
  -isMateMapped true       \
  -isPrimaryAlignment true \
  -in "$tmpdir/$input"     \
  -out "$tmpdir/$output"

if [[ $stopat == 'filter1' ]]; then
  finish "$tmpdir/$output"
fi

echo "--- sort by name started ---"
input="filtered1.bam"
output="sorted2.bam"
java -jar $PICARD_DIR/SortSam.jar \
  I="$tmpdir/$input"              \
  O="$tmpdir/$output"             \
  SO=queryname

# Remove chimeric reads and reads < 100bp
echo "--- dechim and rlen started ---"
input="sorted2.bam"
output="dechim.rlen.bam"  # this script gives a hardcoded output name
python "$scriptdir/rm_chim_in_pair.py" -r $chrom -b $chim_bounds "$tmpdir/$input" "$tmpdir/$output"

if [[ $stopat == 'chimrlen' ]]; then
  finish "$tmpdir/$output"
fi

# bamleftalign
echo "--- realignment started ---"
input="dechim.rlen.bam"
output="realigned.bam"
samtools view -b "$tmpdir/$input" \
  | bamleftalign -f "$ref"        \
  | samtools view -b - > "$tmpdir/$output"

echo "--- sort by coordinate started ---"
input="realigned.bam"
output="sorted3.bam"
java -jar $PICARD_DIR/SortSam.jar \
  I="$tmpdir/$input"              \
  O="$tmpdir/$output"             \
  SO=coordinate
# Also index it.
samtools index "$tmpdir/$output"

if [[ $stopat == 'realign' ]]; then
  finish "$tmpdir/$output"
fi

# Generate consensus sequence
echo "--- major sequence started ---"
input="sorted3.bam"
output="sorted3.bam.major.fa"  # this script gives a hardcoded output name
output2="sorted3.bam.major.fa.fai"
python "$scriptdir/get_major_from_bam.py" "$tmpdir/$input"

# Recalculate NM-tag using new consensus
echo "--- recalculate NM started ---"
input="sorted3.bam"
input2="sorted3.bam.major.fa"
output="nm-recalc.bam"
samtools fillmd -b "$tmpdir/$input" \
  "$tmpdir/$input2"                 \
  2>/dev/null                       \
  > "$tmpdir/$output"

echo "--- sort by name 2 started ---"
input="nm-recalc.bam"
output="sorted4.bam"
java -jar $PICARD_DIR/SortSam.jar \
  I="$tmpdir/$input"              \
  O="$tmpdir/$output"             \
  SO=queryname

# Filter by (recalculated) NM-tag
echo "--- NM selection started ---"
input="sorted4.bam"
output="nm-ratio.sorted4.bam"  # this script gives a hardcoded output name
python "$scriptdir/nm-ratio.select.py" "$tmpdir/$input"

# Re-sort by coordinate, index
echo "--- final sort and index ---"
input="nm-ratio.sorted4.bam"
output="sorted5.bam"
output2="sorted5.bam.bai"
samtools sort "$tmpdir/$input" "$tmpdir/$output.tmp"
mv "$tmpdir/$output.tmp.bam" "$tmpdir/$output"
samtools index "$tmpdir/$output"

finish "$tmpdir/$output"
