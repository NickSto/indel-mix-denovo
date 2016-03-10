#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import argparse
import subprocess

# The locations in the BAM headers that are expected to always be different, even in the correct
# output. Ignore these items.
EXPECTED_DIFFS = {
  '@PG': {  # tag
    1: {      # line
      'CL': {   # key
        'METRICS_FILE', 'OUTPUT', 'INPUT'  # kwargs
      }
    }
  },
  '@00': None  # Use None to indicate that everything under this level should be ignored.
}
ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Compares two BAM files, ignoring differences in the header that don't matter.
The alignments are compared precisely, though, including their order. At the moment it ignores
a few arguments in the CL field of Picard's @PG line. The arguments are paths which can change
without an actual difference in the data."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('expected', metavar='expected.bam',
    help='The expected BAM file.')
  parser.add_argument('actual', metavar='actual',
    help='The actual BAM file.')

  args = parser.parse_args(argv[1:])

  # The files exist?
  if not os.path.isfile(args.expected):
    print('Error: Baseline BAM "{}" not accessible.'.format(args.expected))
  if not os.path.isfile(args.actual):
    print('FAIL: "{}" does not exist or is empty. Pipeline failure?'.format(args.actual))
    return 1

  # Check non-header portion (the reads).
  expected_md5 = get_reads_md5(args.expected)
  actual_md5 = get_reads_md5(args.actual)
  if expected_md5 == actual_md5:
    print('Reads in output BAM are as expected.')
  else:
    print('FAIL: Reads in output BAMs differ!')
    return 1

  # Check header.
  expected_header = get_header(args.expected)
  actual_header = get_header(args.actual)
  diff, diffs = diff_headers(expected_header, actual_header, EXPECTED_DIFFS)
  if diff:
    print('FAIL: Header in output BAM differs in unexpected ways:')
    print(*diffs, sep='\n')
    return 1
  else:
    print('Header in output BAM is exactly as expected (or differs only in the expected ways).')


def get_reads_md5(path):
  try:
    cmd1 = subprocess.Popen(['samtools', 'view', path], stdout=subprocess.PIPE)
    cmd2 = subprocess.Popen(['md5sum'], stdin=cmd1.stdout, stdout=subprocess.PIPE)
    cmd1.stdout.close()
    output = cmd2.communicate()[0]
  except subprocess.CalledProcessError as cpe:
    fail('Command "{}" returned non-zero exit status {}'.format(cpe.cmd, cpe.returncode))
  except OSError:
    fail('Command failed (OSError): $ samtools view "{}" | md5sum'.format(path))
  return output.split()[0]


def get_header(path):
  header = {}
  process = subprocess.Popen(['samtools', 'view', '-H', path], stdout=subprocess.PIPE)
  for line in process.stdout:
    fields = line.strip().split('\t')
    tag = fields[0]
    data = header.get(tag, [])
    datum = {}
    for field in fields[1:]:
      assert field[2] == ':', field
      key = field[:2]
      value = field[3:]
      assert key not in datum, key
      datum[key] = value
    data.append(datum)
    header[tag] = data
  return header


def diff_headers(header1, header2, expected_diffs={}):
  """Find the differences between two headers.
  header1 is the expected, header2 is the actual."""
  diff = False
  diffs = []
  tags = get_intersection(header1.keys(), header2.keys())
  for tag in tags:
    if tag in expected_diffs and expected_diffs[tag] is None:
      pass  # The whole tag is exempted.
    elif tag in header1 and tag in header2:
      # The tag is present in both headers.
      lines1 = header1[tag]
      lines2 = header2[tag]
      i = 0
      while i < len(lines1) or i < len(lines2):
        expected_diffs_lines = expected_diffs.get(tag, {})
        if i in expected_diffs_lines and expected_diffs_lines[i] is None:
          pass  # The whole line is exempted.
        elif i < len(lines1) and i < len(lines2):
          # There is a corresponding line at this position in both headers.
          line1 = lines1[i]
          line2 = lines2[i]
          if format_header_line(tag, line1) == format_header_line(tag, line2):
            pass  # The lines are identical.
          else:
            # The lines differ.
            diffs.append(tag+':')
            # Which key is different?
            keys = get_intersection(line1.keys(), line2.keys())
            for key in keys:
              expected_diffs_keys = expected_diffs_lines.get(i, {})
              if key in expected_diffs_keys and expected_diffs_keys[key] is None:
                pass  # The key is exempted.
              elif line1[key] != line2[key]:
                if key == 'CL':
                  diffs.append('  CL:')
                  nargs1, kwargs1 = parse_command(line1[key])
                  nargs2, kwargs2 = parse_command(line2[key])
                  for narg1, narg2 in zip(nargs1, nargs2):
                    if narg1 != narg2:
                      diff = True
                      diffs.append('  < '+narg1)
                      diffs.append('  > '+narg2)
                  expected_diffs_kwargs = expected_diffs_keys.get(key, {})
                  kwarg_keys = get_intersection(kwargs1.keys(), kwargs2.keys())
                  for kwarg_key in kwarg_keys:
                    if kwarg_key in expected_diffs_kwargs:
                      pass  # The kwarg is exempted.
                    elif kwarg_key in kwargs1 and kwarg_key in kwargs2:
                      if kwargs1[kwarg_key] != kwargs2[kwarg_key]:
                        diff = True
                        diffs.append('  < {}={}'.format(kwarg_key, kwargs1[kwarg_key]))
                        diffs.append('  > {}={}'.format(kwarg_key, kwargs2[kwarg_key]))
                    elif kwarg_key not in kwargs2:
                      diff = True
                      diffs.append('  < {}={}'.format(kwarg_key, kwargs1[kwarg_key]))
                    elif kwarg_key not in kwargs1:
                      diff = True
                      diffs.append('  < {}={}'.format(kwarg_key, kwargs2[kwarg_key]))
                else:
                  diffs.append('< {}:{}'.format(key, line1[key]))
                  diffs.append('> {}:{}'.format(key, line2[key]))
        elif i >= len(lines2):
          # An expected line has disappeared from the output.
          diff = True
          diffs.append('< '+format_header_line(tag, lines1[i]))
        elif i >= len(lines1):
          # An unexpected line has appeared in the output.
          diff = True
          diffs.append('> '+format_header_line(tag, lines2[i]))
        i += 1
    elif tag not in header2:
      # An expected tag is missing from the output.
      for line in header1[tag]:
        diff = True
        diffs.append('< '+format_header_line(tag, line))
    elif tag not in header1:
      # An unexpected tag has appeared in the output.
      for line in header2[tag]:
        diff = True
        diffs.append('> '+format_header_line(tag, line))
  return diff, diffs


def get_intersection(l1, l2):
  """Get the intersection of two lists (as a list)."""
  s1 = set(l1)
  s2 = set(l2)
  set_intersection = s1.intersection(s2)
  return list(set_intersection)


def format_header_line(tag, line):
  keys = sorted(line.keys())
  fields = [key+':'+line[key] for key in keys]
  return tag+'\t'+'\t'.join(fields)


def parse_command(command):
  args = command.split()
  positionals = []
  kwargs = {}
  for arg in args:
    # Is it a Picard/GATK-style keyword argument?
    match = re.search(r'^([A-Z0-9_]+)=(.*)$', arg)
    if match:
      kwargs[match.group(1)] = match.group(2)
    else:
      positionals.append(arg)
  return positionals, kwargs


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
