#!/usr/bin/env python
# coding: utf-8

import os
import sys
seq2taxid_file = '/mnt/c/Weekly_Miseq/rRNA/seq2fulltax.tsv'
no_hit_file = "/mnt/c/Weekly_MiSeq/rRNA/No_hits.txt"

prev_seqs = {}
# get prematched taxids
for file in (seq2taxid_file, no_hit_file):
    if os.path.isfile(file):
        with open(file, "r") as in_fh:
            for line in in_fh:
                line = line.strip().split("\t")
                prev_seqs[line[0]] = 1

# get sequences
sequences = []
for file in os.listdir():
    if file.endswith("uc.fa"):
        with open(file, "r") as in_fh:
            seqs = in_fh.read().split(">")
            sequences += ["".join(seq.split("\n")[1:]) for seq in seqs if seq]

print(len(sequences))
sequences = set(sequences)
print(len(sequences))

# make file for blastn
cs = 0
with open("seqs4blast.fa", "w") as out_fh:
    counter = 0
    for seq in sequences:
        if not seq in prev_seqs:
            out_fh.write(f">{counter}\n{seq}\n")
            counter += 1
print(counter)
