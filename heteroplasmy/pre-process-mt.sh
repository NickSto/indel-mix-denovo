#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

# note: the reference does not have to be indexed; bamleftalign will do that.
REF_DEFAULT="$HOME/bx/data/chrM-rCRS.fa"
CHROM_DEFAULT="chrM"
CHROM_LEN_DEFAULT=""
MARGIN_DEFAULT=""
BOUNDS_DEFAULT="600 16000"
MIN_LEN_DEFAULT="100"
REQUIRED_SCRIPTS="rm_chim_in_pair.py nm-ratio.select.py get_major_from_bam.py"
REQUIRED_COMMANDS="dirname basename file java bamtools samtools bamleftalign"
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
-c: A chromosome to target the analysis to. Default: $CHROM_DEFAULT. If you change this, you'll
    probably need to provide either -B or -m since the default -B ($BOUNDS_DEFAULT) will probably
    no longer be valid.
    N.B.: Not fully implemented yet! Only works up through chimrlen!
-m: The length of the ends of the chromosome that should be exempt from the requirement that reads
    be non-chimeric. This is to allow for circular chromosomes, so that reads that span the start
    and end aren't mistakenly discarded as chimeric.
-B: Bounds to hand to rm_chim_in_pair.py. Default: $BOUNDS_DEFAULT
-L: The length of the target chromosome. Should not be necessary, as the script should be able to
    determine this using the reference file and the target chromosome name (-c).
-M: The minimum read length to accept when filtering reads in rm_chim_in_pair.py. Default: $MIN_LEN_DEFAULT"

function fail {
  echo "Error: $1" >&2
  exit 1
}

# Create an equivalent of readlink -f for systems that lack it (BSD).
function readlink_f {
  if readlink -f dummy >/dev/null 2>/dev/null; then
    readlink -f "$1"
  else
    perl -MCwd -le 'print Cwd::abs_path(shift)' "$1"
  fi
}

# Get this script's actual directory, resolving links.
scriptdir=$(dirname $(readlink_f ${BASH_SOURCE[0]}))

# read in arguments
ref="$REF_DEFAULT"
chrom="$CHROM_DEFAULT"
margin="$MARGIN_DEFAULT"
chrom_len="$CHROM_LEN_DEFAULT"
chim_bounds="$BOUNDS_DEFAULT"
min_len="$MIN_LEN_DEFAULT"
stopat=''
while getopts ":r:s:c:B:m:L:M:h" opt; do
  case "$opt" in
    r) ref="$OPTARG";;
    s) stopat="$OPTARG";;
    c) chrom="$OPTARG";;
    B) chim_bounds="$OPTARG";;
    m) margin="$OPTARG";;
    L) chrom_len="$OPTARG";;
    M) min_len="$OPTARG";;
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
# If -L was not given, need to use either bioawk or seqlen.awk to determine target sequence length.
if ! [[ $chrom_len ]]; then
  if ! which bioawk >/dev/null 2>/dev/null; then
    REQUIRED_SCRIPTS="$REQUIRED_SCRIPTS seqlen.awk"
  fi
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
type=$(file -b --mime-type "$(readlink_f "$inpath")")
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

# If -L not was given, use bioawk or seqlen.awk to determine target chromosome length.
if ! [[ $chrom_len ]]; then
  if which bioawk >/dev/null 2>/dev/null; then
    chrom_len=$(bioawk -c fastx '$name == "'$chrom'" { print length($seq) }' $ref | head -n 1)
  else
    chrom_len=$(awk -f $scriptdir/seqlen.awk -v target=$chrom $ref | cut -f 2 | head -n 1)
  fi
fi

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
  VALIDATION_STRINGENCY=SILENT    \
  I="$input"                      \
  O="$tmpdir/$output"             \
  SO=coordinate

# Mark duplicates
echo "--- markDup started ---"
input="sorted1.bam"
output="dupmarked.bam"
output2="dupmarked.metrics"
java -jar $PICARD_DIR/MarkDuplicates.jar \
  VALIDATION_STRINGENCY=SILENT           \
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
  VALIDATION_STRINGENCY=SILENT    \
  I="$tmpdir/$input"              \
  O="$tmpdir/$output"             \
  SO=queryname

# Remove chimeric reads and reads < $min_len bp
#TODO: Instead of giving $min_len directly with -M, allow giving an input read
#      length and minimum length in percent.
echo "--- dechim and rlen started ---"
input="sorted2.bam"
output="dechim.rlen.bam"
if [[ $chrom_len ]] && [[ $margin ]]; then
  python "$scriptdir/rm_chim_in_pair.py" \
    -r $chrom                            \
    -l $min_len                          \
    -m $margin                           \
    -L $chrom_len                        \
    "$tmpdir/$input"                     \
    "$tmpdir/$output"
else
  python "$scriptdir/rm_chim_in_pair.py" \
    -r $chrom                            \
    -l $min_len                          \
    -b $chim_bounds                      \
    "$tmpdir/$input"                     \
    "$tmpdir/$output"
fi

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
  VALIDATION_STRINGENCY=SILENT    \
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
output="sorted3.bam.major.fa"
python "$scriptdir/get_major_from_bam.py" \
  -c $chrom                               \
  -l $chrom_len                           \
  "$tmpdir/$input"                        \
  -o "$tmpdir/$output"

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
  VALIDATION_STRINGENCY=SILENT    \
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
