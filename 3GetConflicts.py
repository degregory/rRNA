#!/usr/bin/env python
# coding: utf-8

import os
import sys

taxid2tax_file = '/mnt/c/Weekly_MiSeq/rRNA/taxid2fulltax.tsv'
all12s = "all12S.blastn.nt.tsv"
seq2taxid_file = '/mnt/c/Weekly_Miseq/rRNA/seq2fulltax.tsv'
no_hit_file = "/mnt/c/Weekly_MiSeq/rRNA/No_hits.txt"
resolved_file = '/mnt/c/Weekly_MiSeq/rRNA/conflicts.resolved.tsv'



prev_seqs = {}
# get prematched taxids
for file in (seq2taxid_file, no_hit_file):
    if os.path.isfile(file):
        with open(file, "r") as in_fh:
            for line in in_fh:
                line = line.strip().split("\t")
                prev_seqs[line[0]] = 1

blast_dict = {}
samp_dict = {}

for file in os.listdir():
    if file.endswith(".uc.fa"):
        samp =  file.split(".")[0]
        with open(file, "r") as in_fh:
            seqs = in_fh.read().split(">")
            for seq in seqs:
                if seq:
                    seq = seq.split("\n")
                    seq[1] = "".join(seq[1:])
                    try:
                        samp_dict[seq[1]].append(samp)
                    except KeyError:
                        samp_dict[seq[1]] = [samp]
print(len(samp_dict))

blastee_file = "seqs4blast.fa"
seq_dict = {}
with open(blastee_file, "r") as in_fh:
    seqs = in_fh.read().split(">")
    for seq in seqs:
        if seq:
            seq = seq.split("\n")
            seq[1] = "".join(seq[1:])
            seq_dict[seq[0]] = seq[1]
print(len(seq_dict))

blasted_file = 'blastn.tsv'
with open(blasted_file, "r") as in_fh:
    for line in in_fh:
        line = line.split("\t")
        if line[0] in seq_dict:
            try:
                blast_dict[seq_dict[line[0]]].append("\t".join(line[1:]))
            except KeyError:
                blast_dict[seq_dict[line[0]]] = ["\t".join(line[1:])]

print(len(blast_dict))

# get local taxid and associated full tax
fulltax_dict = {}
if os.path.isfile(taxid2tax_file):
    with open(taxid2tax_file, "r") as in_fh:
        for line in in_fh:
            line = line.split("\t")
            fulltax_dict[line[0]] = line[1].strip()

prev_res = {}
if os.path.isfile(resolved_file):
    with open(resolved_file, "r") as in_fh:
        for line in in_fh:
            prev_res[line.split("\t")[1]] = 1

with open(all12s, "a") as out_fh, open("conflicts.tsv", "w") as conflicts:
    conflicts.write("\tsequence\tsamples\tstaxid\tsskingdom\tsscinames\tscomname\tcount\n")
    counter = 0
    for seq, bl_lines in blast_dict.items():
        if not seq in prev_seqs:
            for line in bl_lines:
                out_fh.write(f"{seq}\t{line}")
            taxes = [bl_info.split("\t")[8] for bl_info in bl_lines]
            if not len(set(taxes)) == 1 and not seq in prev_res:
                if "Eukaryota" in [bl_info.split("\t")[10] for bl_info in bl_lines]:
                    p = False
                    checked_lines = []
                    for line in bl_lines:
                        if "\t100.000\t0\t100\t" in line and int(line.split("\t")[7]) == len(seq):
                            checked_lines.append(line)
                    if len(checked_lines) > 1:
                        if "Eukaryota" in [bl_info.split("\t")[10] for bl_info in checked_lines]:
                            taxinfs = {}
                            domains = {}
                            for line in checked_lines:
                                line = line.split("\t")
                                tax_info = f"{line[8]}\t{line[10]}\t{line[11]}\t{line[13]}"
                                try:
                                    taxinfs[tax_info] += 1
                                except KeyError:
                                    taxinfs[tax_info] = 1
                                try:
                                    domains[line[10]] += 1
                                except KeyError:
                                    domains[line[10]] = 1
                            if len(taxinfs) > 1 and domains["Eukaryota"] > .05 * sum(domains.values()):
                                
                                for taxinf in taxinfs:
                                    conflicts.write(f"{counter}\t{seq}\t{' '.join(samp_dict[seq])}\t{taxinf}\t{taxinfs[taxinf]}\t")
                                    try:
                                        conflicts.write(fulltax_dict[taxinf.split("\t")[0]])
                                    except:
                                        pass
                                    conflicts.write("\n")
                                counter += 1
    print(counter)

nh_count = 0
with open(no_hit_file, "a") as no_hits:
    for seq in samp_dict:
        if not (seq in blast_dict or seq in prev_seqs):
            nh_count += 1
            no_hits.write(f"{seq}\n")
print(nh_count)


seqs = []
with open(no_hit_file, "r") as no_hits:
    seqs = set(no_hits.read().split("\n"))


with open(no_hit_file, "w") as no_hits:
    for seq in seqs:
        no_hits.write(f"{seq}\n")
