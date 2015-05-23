#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

USAGE="\$ $(basename $0) families.tsv [pipeline-family.py args]
This is a wrapper script to launch multiple instances of pipeline-family.py in
the background, in parallel.
Provide the same arguments you would for pipeline-family.py, but without
-s and -f, and inserting at the front the filename of the families table. This
should be a space-delimited file with one family per line. Each line should have
the family id as the first field, and the rest of the fields should be the ids
of the samples in that family."


function main {

  if [[ $# -eq 0 ]] || [[ $1 == '-h' ]]; then
    fail "$USAGE"
  fi

  # Put $1 in $families and shift it off the front of $@ so we can use the rest as the arguments to
  # pipeline-family.py.
  families="$1"
  shift

  script_dir=$(real_dir)
  if ! [[ -f "$families" ]]; then
    fail "Error: families table $families.tsv missing."
  fi
  if ! [[ -f "$script_dir/pipeline-family.py" ]]; then
    fail "Error: missing script: $script_dir/pipeline-family.py"
  fi

  while read family samples; do
    if [[ ! "$family" ]] || [[ ! "$samples" ]]; then
      continue
    fi
    sample_list=$(echo $samples | tr ' ' ',')
    echo '+ $' python "$script_dir/pipeline-family.py" -f $family -s $sample_list "$@"
    python "$script_dir/pipeline-family.py" -f $family -s $sample_list "$@" &
  done < "$families"

}


# Get the script's actual directory path
function real_dir {
  # Does readlink -f work? (It doesn't on BSD.)
  if readlink -f dummy >/dev/null 2>/dev/null; then
    dirname $(readlink -f ${BASH_SOURCE[0]})
  else
    # If readlink -f doesn't work (like on BSD).
    # Read the link destination from the output of ls -l and cd to it.
    # Have to cd to the link's directory first, to handle relative links.
    # With help from https://stackoverflow.com/a/246128/726773
    unset CDPATH
    local source="${BASH_SOURCE[0]}"
    while [[ -h "$source" ]]; do
      local dir="$(cd -P $(dirname "$source") && pwd)"
      local link="$(ls -l "$source" | awk '{print $NF}')"
      # absolute or relative path?
      if [[ "$link" == /* ]]; then
        source="$link"
      else
        source="$dir/$link"
      fi
    done
    dir="$(cd -P $(dirname "$source") && pwd)"
    echo "$dir"
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
