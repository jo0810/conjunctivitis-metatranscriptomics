"""Shared loading and filtering of read-level ARG alignment metrics."""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]          # this folder and the repository root
from config import in_file  # noqa: E402

BREADTH_MIN = 50    # coverage breadth: % of reference protein residues covered by >= 1 read
DPM_MIN = 1         # an ARG is retained in a library if DPM > 1 ...
READS_MIN = 10      # ... or read count > 10


def libraries():
    """Analysed libraries, in the order listed in library_metadata.csv."""
    meta = pd.read_csv(in_file('library_metadata.csv'))
    return [str(s) for s in meta['sample']], meta


def filtered_alignments(libs):
    """Read-level ARG metrics restricted to the analysed libraries and to retained ARGs."""
    a = pd.read_csv(in_file('arg_read_level_metrics.tsv'), sep='\t')
    a = a[a['sample'].astype(str).isin(libs)].copy()
    p = a[(a['coverage_breadth_pct'] >= BREADTH_MIN) &
          ((a['depth_per_million'] > DPM_MIN) | (a['read_count'] > READS_MIN))].copy()
    return p


def annotation():
    """CARD-based annotation of ARG features, indexed by feature_id."""
    return pd.read_csv(in_file('arg_feature_annotation.tsv'), sep='\t').set_index('feature_id')
