#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)


function main {
  if [[ $# -eq 0 ]] || [[ $1 == all ]]; then
    pipeline
    return
  fi

  # Run the requested tests
  for test in "$@"; do
    case "$test" in
      pipeline)
        pipeline;;
      pipeline2)
        pipeline2;;
      pipeline3)
        pipeline3;;
    esac
  done
}


########## Functional tests ##########

# Do all tests of the pipeline.
function pipeline {
  pipeline2
  pipeline3
}

# Test up through step 2 of the pipeline.
function pipeline2 {
  echo -e "\tpipeline.py steps 1 and 2 ::: R39-M249-reduced.bam:"
  mkdir "$dirname/pipeline/tmp" || return
  python "$dirname/../pipeline.py" -E 2 -s M249 -R chrM "$dirname/chrM-rCRS.fa" \
    "$dirname/pipeline/R39-M249-reduced_1.fq" "$dirname/pipeline/R39-M249-reduced_2.fq" \
    "$dirname/pipeline/tmp" 2>/dev/null
  #TODO: edit diff.sh to allow different paths in @PG command line.
  bash "$dirname/diff.sh" "$dirname/pipeline/out2.R39-M249-reduced.bam" \
    "$dirname/pipeline/tmp/bam1filt.bam"
  rm -r "$dirname/pipeline/tmp"
}

# Test step 3 of the pipeline.
function pipeline3 {
  echo -e "\tpipeline.py step 3 ::: R39-M249-reduced.bam:"
  mkdir "$dirname/pipeline/tmp" || return
  cp "$dirname/pipeline/out2.R39-M249-reduced.bam" "$dirname/pipeline/tmp/bam1filt.bam"
  python "$dirname/../pipeline.py" -B 3 -E 3 -s M249 -R chrM "$dirname/chrM-rCRS.fa" \
    "$dirname/pipeline/R39-M249-reduced_1.fq" "$dirname/pipeline/R39-M249-reduced_2.fq" \
    "$dirname/pipeline/tmp" 2>/dev/null
  #TODO: edit diff.sh to allow different paths in @PG command line.
  bash "$dirname/diff.sh" "$dirname/pipeline/out3.R39-M249-reduced.bam" \
    "$dirname/pipeline/tmp/bam1dedup.bam"
  rm -r "$dirname/pipeline/tmp"
}

main "$@"
