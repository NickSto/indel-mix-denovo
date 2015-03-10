#!/usr/bin/env bash

function main {
  if [[ $# -lt 2 ]] || [[ $1 == '-h' ]]; then
    echo "Usage: \$ $(basename $0) expected.bam actual.bam
Compares two BAM files, ignoring differences in the header that don't matter.
The alignments are compared precisely, though, including their order. At the
moment it just ignores the @PG lines in the header, which include the command
line, which can change depending on the pwd at the time of invocation." >&2
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
    return 1
  fi

  # check non-header portion (the reads)
  expected_md5=$(md5sum <(samtools view "$expected") | awk '{print $1}')
  actual_md5=$(md5sum <(samtools view "$actual") | awk '{print $1}')
  if [[ $expected_md5 == $actual_md5 ]]; then
    echo "Reads in output BAM are as expected."
  else
    echo "FAIL: Reads in output BAMs differ!"
    return 1
  fi

  # check header
  diff=$(samtools view -H "$expected" | diff - <(samtools view -H "$actual"))
  if [[ "$diff" ]]; then
    diff_tags=$(echo "$diff" | awk '$1 == ">" || $1 == "<" {print $2}')
    # $expected should be false if any of the differing tags isn't an expected one.
    expected=true
    for tag in $diff_tags; do
      # NOTE: Add header tags here!
      case $tag in
        @PG)
          true;;
        *)
          expected='';;
      esac
    done
    if [[ $expected ]]; then
      echo "Header in output BAM differs in only the expected ways."
    else
      echo "FAIL: Header in output BAM differs in unexpected ways:"
      echo "$diff"
      return 1
    fi
  else
    echo "Header in output BAM is identical to the expected one."
  fi
}

main "$@"
