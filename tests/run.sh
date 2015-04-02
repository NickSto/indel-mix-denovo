#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)
tmp="$dirname/pipeline/test.tmp"

USAGE="Usage: \$ $(basename $0) [options] [test1 [test2]]"


function main {

  do_all=true
  verbose=''
  # Run the requested tests
  for arg in "$@"; do
    # Check for options
    #TODO: option to keep test data at end instead of removing it.
    if [[ ${arg:0:1} == '-' ]]; then
      case "$arg" in
        -h)
          echo "$USAGE" >&2
          echo "Currently valid tests:" >&2
          list_tests >&2
          exit 1;;
        -v)
          verbose=true;;
        *)
          echo "Unrecognized option \"$arg\"." >&2;;
      esac
      continue
    fi
    # Execute valid tests (if they're existing functions).
    if [[ $(type -t $arg) == function ]]; then
      do_all=''
      if [[ $verbose ]]; then
        $arg
      else
        $arg 2>/dev/null
      fi
    else
      echo "Unrecognized test \"$arg\"." >&2
      do_all=''
    fi
  done

  # If no tests were specified in arguments, do all tests.
  if [[ $do_all ]]; then
    if [[ $verbose ]]; then
      all
    else
      all 2>/dev/null
    fi
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

function failout {
  echo "$@"
  exit 1
}

function list_tests {
  while read declare f test; do
    # Filter out functions that aren't tests.
    if echo "$initial_declarations" | grep -qF 'declare -f '"$test"; then
      continue
    else
      echo "$test"
    fi
  done < <(declare -F)
}

# Capture a list of all functions defined before the tests, to tell which are actual functions
# and which are tests.
initial_declarations=$(declare -F)

########## Functional tests ##########

# Do all tests.
function all {
  pipeline
}

# Do all pipeline tests.
function pipeline {
  pipeline2 && echo
  pipeline3 && echo
  pipefull11 && echo
  pipefull
}

# Test up through step 2 of the pipeline.
function pipeline2 {
  echo -e "\tpipeline.py steps 1 and 2 ::: R39-M249-reduced_[12].fq:"
  mkdir "$tmp" || failout "Error: tmp dir $tmp exists."
  python "$dirname/../pipeline.py" -E 2 -s M249 -r chrM -l 250 "$dirname/chrM-rCRS.fa" \
    "$dirname/pipeline/R39-M249-reduced_1.fq" "$dirname/pipeline/R39-M249-reduced_2.fq" "$tmp"
  #TODO: edit diff.sh to allow different paths in @PG command line.
  bash "$dirname/bam-diff.sh" "$dirname/pipeline/out2.R39-M249-reduced.bam" \
    "$tmp/bam1filt.bam"
  rm -r "$tmp"
}

# Test step 3 of the pipeline.
function pipeline3 {
  echo -e "\tpipeline.py step 3 ::: R39-M249-reduced_[12].fq:"
  mkdir "$tmp" || failout "Error: tmp dir $tmp exists."
  cp "$dirname/pipeline/out2.R39-M249-reduced.bam" "$tmp/bam1filt.bam"
  python "$dirname/../pipeline.py" -B 3 -E 3 -s M249 -r chrM -l 250 "$dirname/chrM-rCRS.fa" \
    "$dirname/pipeline/R39-M249-reduced_1.fq" "$dirname/pipeline/R39-M249-reduced_2.fq" "$tmp"
  #TODO: edit diff.sh to allow different paths in @PG command line.
  bash "$dirname/bam-diff.sh" "$dirname/pipeline/out3.R39-M249-reduced.bam" \
    "$tmp/bam1dedup.bam"
  rm -r "$tmp"
}

# Test the pipeline from steps 1-11
function pipefull11 {
  echo -e "\tpipeline.py steps 1-11 ::: R39-M249-reduced_[12].fq:"
  mkdir "$tmp" || failout "Error: tmp dir $tmp exists."
  python "$dirname/../pipeline.py" -E 11 -s M249 -r chrM -l 250 "$dirname/chrM-rCRS.fa" \
    "$dirname/pipeline/R39-M249-reduced_1.fq" "$dirname/pipeline/R39-M249-reduced_2.fq" "$tmp"
  # Diff the results. The date ends up in the final vcf, so omit that.
  exceptions='^##(reference|fileDate)='
  grep -Ev "$exceptions" "$dirname/pipeline/out11.R39-M249-reduced.vcf" \
    | diff - <(grep -Ev "$exceptions" "$tmp/nvc.vcf") > "$tmp/diff.txt"
  lines=$(wc -l "$tmp/diff.txt" | awk '{print $1}')
  if [[ $lines == 0 ]]; then
    echo "Output nvc.vcf and out11.R39-M249-reduced.vcf differ only by the expected amount."
  elif ! [[ -e "$tmp/nvc.vcf" ]]; then
    echo "FAIL: Output file nvc.vcf missing."
  elif [[ $lines -le 40 ]]; then
    cat "$tmp/diff.txt"
    echo "FAIL: Output nvc.vcf and expected out11.R39-M249-reduced.vcf differ."
  else
    echo "FAIL: Output nvc.vcf and expected out11.R39-M249-reduced.vcf differ. Diff lines: $lines"
  fi
  rm -r "$tmp"
}

# Test the pipeline all the way through
function pipefull {
  echo -e "\tpipeline.py all steps ::: G3825.1a_[12].fq:"
  mkdir "$tmp" || failout "Error: tmp dir $tmp exists."
  python "$dirname/../pipeline.py" -s G3825.1a -r EBOV -l 100 -c 100 -m 0 "$dirname/EBOV.fa" \
    "$dirname/pipeline/G3825.1a_1.fq" "$dirname/pipeline/G3825.1a_2.fq" "$tmp"
  diff -s "$dirname/pipeline/out14.G3825.1a.tsv" "$tmp/vars.tsv"
  rm -r "$tmp"
}

main "$@"
