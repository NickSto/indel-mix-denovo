BEGIN {
  FS="\t"
  OFS="\t"
  # Command line options (set with -v)
  if (! covg) {
    covg = 1000
  }
  if (! freq) {
    freq = 1
  }
  if (! strand) {
    strand = 1.0
  }
  if (! mate) {
    mate = 1.0
  }
  # Format for excluded: regions separated by commas, start/end coordinates separated by dashes.
  # Chromosome names can be given with coordinates, joined by a colon. But chromosome names must not
  # contain a dash.
  # e.g. "302-310,16183-16192" or "chrM:302-310" or "chrM:302-chrM:310,16183-16192"
  if (excluded) {
    # Split the comma-separated regions.
    split(excluded, excluded_strs, ",")
    i = 0
    for (j in excluded_strs) {
      # Split the start-end coordinates.
      split(excluded_strs[j], region, "-")
      # Split the starting chromosome name, if any, from the coordinate.
      split(region[1], fields, ":")
      if (fields[2]) {
        start_chrs[i]   = fields[1]
        start_coords[i] = fields[2]
      } else {
        start_chrs[i]   = ""
        start_coords[i] = fields[1]
      }
      # Split the ending chromosome name, if any, from the coordinate.
      split(region[2], fields, ":")
      if (fields[2]) {
        end_chrs[i]   = fields[1]
        end_coords[i] = fields[2]
      } else {
        end_chrs[i]   = ""
        end_coords[i] = fields[1]
      }
      i++
    }
  }
}

$6 > covg && $8 > freq && $20 < strand && $21 < mate {
  # Filter out excluded regions, too.
  passed = 1
  for (i in start_coords) {
    if (start_chrs[i]) {
      if (start_chrs[i] == $2 && $3 >= start_coords[i] && $3 <= end_coords[i]) {
        passed = 0
      }
    } else {
      if ($3 >= start_coords[i] && $3 <= end_coords[i]) {
        passed = 0
      }
    }
  }
  if (passed) {
    print $0
  }
}
