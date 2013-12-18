#!/usr/bin/env python
"""Methods to select reads from a BAM that contain given variants, and return
statistics on their characteristics."""
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
import sys