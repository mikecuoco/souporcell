#!/usr/bin/env python
# Author: Mike Cuoco
# Created on: 6/24/22, 5:21 PM
#
# Description: Split clusters.tsv from generated from merged bam into separate files corresponding to original bams
# Usage: python {filename} <arg1> <arg2>

import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description="Split clusters.tsv from generated from merged bam into separate files corresponding to original bams")
parser.add_argument(
    "--clusters", type=Path, help="path to clusters.tsv file from souporcell"
)
parser.add_argument(
    "--merge-info", type=Path, help="path to 3-column {name}_info.tsv, compatible output from merge_10x_bams.py"
)

args = parser.parse_args()

info = pd.read_csv(args.merge_info, sep="\t").set_index("CB_pool")
clusters = pd.read_csv(args.clusters, sep="\t").set_index("barcode")
merged = pd.merge(info,clusters, how="outer", left_index=True, right_index=True)

for id,run in merged.groupby('Run_id'):

    run.drop(labels = 'Run_id', axis = 1, inplace = True)
    run.rename(columns = {"CB_origin":"barcode"}, inplace = True)

    # write out clusters for this run
    run.to_csv(args.clusters.parent / (id + "_clusters.tsv"), sep="\t", index=False)
    