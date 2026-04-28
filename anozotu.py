#!/bin/env python3

import os
import sys
import argparse


parser = argparse.ArgumentParser(
        description=' '
)

parser.add_argument(
    'samp',
    type=str,
    help='prefix'
)

args = parser.parse_args()

samp = args.samp
count_dict = {}
with open(f"{samp}.unoise3.txt", "r") as in_fh:
    for line in in_fh:
        line = line.strip().split("\t")
        if line[1] == "denoise" and line[2].startswith("amp"):
            count_dict[line[2][3:]] = line[0].split("=")[-1]

with open(f"{samp}.un.fa", "r") as in_fh, open(f"{samp}.uc.fa", "w") as out_fh:
    for line in in_fh:
        if line.startswith(">"):
            line = line.strip()
            line += f";size={count_dict[line.split('u')[-1]]}\n"
        out_fh.write(line)