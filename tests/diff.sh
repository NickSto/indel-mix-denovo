#!/usr/bin/env bash

function main {
  if [[ $# -lt 2 ]]; then
    echo "Usage: \$ $(basename $0) expected.bam actual.bam
Compares two BAMs produced by the pre-process-mt.sh filtering pipeline.
Determines if they differ beyond the expected differences.
The expected difference this allows for is the BAM header Picard inserts, which
includes the input filename (which will change for different runs)." >&2
    exit 1
  else
    pipeline "$1" "$2"
  fi
}


# Arg1 is the path to the expected file.
# Arg2 is the path to the file in question.
function pipeline {
  expected="$1"
  actual="$2"

  # the files exist?
  if [[ ! -s "$expected" ]]; then
    echo "Error: Baseline BAM \"$expected\" does not exist or is empty."
    return
  fi
  if [[ ! -s "$actual" ]]; then
    echo "FAIL: \"$actual\" does not exist or is empty. Pipeline failure?"
    return
  fi

  # check non-header portion (the reads)
  expected_crc=$(crc32 <(samtools view "$expected"))
  actual_crc=$(crc32 <(samtools view "$actual"))
  if [[ $expected_crc == $actual_crc ]]; then
    echo "Reads in output BAM are as expected."
  else
    echo "FAIL: Reads in output BAMs differ!"
    return
  fi

  # check header
  diff1=$(samtools view -H "$expected" | diff - <(samtools view -H "$actual"))
  if [[ $diff1 ]] && [[ $(echo "$diff1" | wc -l) == 4 ]]; then
    # get a sub-diff: which fields within the line differ
    line1=$(echo "$diff1" | head -n 2 | tail -n 1 | tr -s ' \t' '\n\n')
    line2=$(echo "$diff1" | head -n 4 | tail -n 1 | tr -s ' \t' '\n\n')
    diff2=$(echo "$line1" | diff - <(echo "$line2"))
    if [[ $(echo "$line1" | head -n 2 | tail -n 1) != "@PG" ]]; then
      echo "FAIL: Header differs only in one line, but not the correct one:"
      echo "$diff1"
    elif [[ $(echo "$diff2" | wc -l) == 12 ]]; then
      inputs=$(echo "$diff2" | grep -c -E '^(<|>) INPUT=')
      outputs=$(echo "$diff2" | grep -c -E '^(<|>) OUTPUT=')
      metrics=$(echo "$diff2" | grep -c -E '^(<|>) METRICS_FILE=')
      if [[ $inputs == 2 && $outputs == 2 && $metrics == 2 ]]; then
        # success! the diff looks as expected!
        echo "Header in output BAM is as expected."
      else
        echo "FAIL: Header differs only in the correct line, but not in the expected way:"
        echo "$diff2"
        return
      fi
    else
      echo "FAIL: Header differs only in the correct line, but not in the expected way:"
      echo "$diff2"
      return
    fi
  elif [[ ! $diff1 ]]; then
    # no difference is good!
    echo "Header in output BAM is as expected."
  else
    # there is a difference, and it's the expected one.
    echo "FAIL: Header in output BAM not as expected:"
    echo "$diff1"
    return
  fi
}

main "$@"
