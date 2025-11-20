#!/usr/bin/env python
# coding: utf-8

import os
import sys

seq2taxid_file = 'all12S.blastn.nt.tsv'
resolved_file = '/mnt/c/Weekly_MiSeq/rRNA/conflicts.resolved.tsv'
taxid2tax_file = '/mnt/c/Weekly_MiSeq/rRNA/taxid2fulltax.tsv'
seq2fulltax_file = '/mnt/c/Weekly_MiSeq/rRNA/seq2fulltax.tsv'


# get existing full tax
seq2fulltax_dict = {}
if os.path.isfile(seq2fulltax_file):
    with open(seq2fulltax_file, "r") as in_fh:
        for line in in_fh:
            line = line.strip().split("\t")
            seq2fulltax_dict[line[0]] = (line[1], line[2], line[3])
print(len(seq2fulltax_dict))

# get prematched taxids
seq2tax_dict = {}
if os.path.isfile(seq2taxid_file):
    with open(seq2taxid_file, "r") as in_fh:
        cur_seq = ["", {}]
        for line in in_fh:
            sline = line.strip().split("\t")
            if sline[0] in seq2fulltax_dict:
                continue
            elif cur_seq[0] and cur_seq[0] == sline[0]:
                if "p" in cur_seq[1]:
                    if "\t100.000\t0\t100\t" in line and int(sline[8]) == len(sline[0]):
                        cur_seq[1]["p"].append(sline[9])
                else:
                    if "\t100.000\t0\t100\t" in line and int(sline[8]) == len(sline[0]):
                        cur_seq[1]["p"] = [sline[9]]
                    else:
                        try:
                            cur_seq[1][int(sline[8])].append(sline[9])
                        except KeyError:
                            cur_seq[1][int(sline[8])] = [sline[9]]
            else:
                if cur_seq[0]:
                    seq2tax_dict[cur_seq[0]] = [] ### track perfect hits
                    if 'p' in cur_seq[1]:
                        seq2tax_dict[cur_seq[0]] = ["p", cur_seq[1]["p"]]
                    else:
                        scores = sorted(cur_seq[1].keys(), reverse=True)[:2]
                        seq2tax_dict[cur_seq[0]] += cur_seq[1][scores[0]]
                        if len(scores) > 1 and scores[1] > scores[0] - 10:
                            seq2tax_dict[cur_seq[0]] += cur_seq[1][scores[1]]
                        # seq2tax_dict[cur_seq[0]] = set(seq2tax_dict[cur_seq[0]])
                cur_seq = [sline[0], {}]
                if "\t100.000\t0\t100\t" in line and int(sline[8]) == len(sline[0]):
                    cur_seq[1]["p"] = [sline[9]]
                else:
                    cur_seq[1][int(sline[8])] = [sline[9]]

        if cur_seq[0]:
            seq2tax_dict[cur_seq[0]] = []
            if 'p' in cur_seq[1]:
                seq2tax_dict[cur_seq[0]] = ["p", cur_seq[1]["p"]]
            else:
                scores = sorted(cur_seq[1].keys(), reverse=True)[:2]
                seq2tax_dict[cur_seq[0]] += cur_seq[1][scores[0]]
                if len(scores) > 1 and scores[1] > scores[0] - 10:
                    seq2tax_dict[cur_seq[0]] += cur_seq[1][scores[1]]

if os.path.isfile(resolved_file):
    with open(resolved_file, "r") as in_fh:
        in_fh.readline()
        cur_seq = ["", []]
        for line in in_fh:
            sline = line.strip().split("\t")
            if sline[1] in seq2fulltax_dict:
                continue
            elif cur_seq[0] and cur_seq[0] == sline[1]:
                cur_seq[1].append(sline[3])
            else:
                if cur_seq[0]:
                    seq2tax_dict[cur_seq[0]] = ["p", cur_seq[1]]
                cur_seq = [sline[1], [sline[3]]]

        if cur_seq[0]:
            seq2tax_dict[cur_seq[0]] = ["p", cur_seq[1]]

print(len(seq2tax_dict))

# get local taxid and associated full tax
fulltax_dict = {}
if os.path.isfile(taxid2tax_file):
    with open(taxid2tax_file, "r") as in_fh:
        for line in in_fh:
            line = line.split("\t")
            if "domain_Bacteria" in line[1]:
                fulltax_dict[line[0]] = "Bacteria"
            elif "domain_Archaea" in line[1]:
                fulltax_dict[line[0]] = "Archaea"
            elif line[1].startswith("no rank_"):
                fulltax_dict[line[0]] = "No Rank"
            else:
                fulltax_dict[line[0]] = line[1].strip()

if seq2tax_dict:
    for seq, taxids in seq2tax_dict.items():
        taxes = {}
        perf = ''
        mixed = ''
        if taxids[0] == 'p':
            perf = 'p'
            taxids = taxids[1]
        for taxid in taxids:
            try:
                fulltax_dict[taxid]
                try:
                    taxes[fulltax_dict[taxid]] += 1
                except:
                    taxes[fulltax_dict[taxid]] = 1
            except KeyError:
                if not taxid == "N/A":
                    print(f"{taxid} not found")
                try:
                    taxes[taxid] += 1
                except KeyError:
                    taxes[taxid] = 1
        if len(taxes) == 1:
            taxes = list(taxes.keys())[0]
            if taxes.startswith("cellular root_cellular organisms;"):
                taxes = "; ".join(taxes.split("; ")[1:])
        else:
            if "N/A" in taxes and len(taxes) > 1:
                del taxes["N/A"]
            if "No Rank" in taxes and len(taxes) > 1:
                del taxes["No Rank"]
            if "Bacteria" in taxes and taxes["Bacteria"] < .05 * sum(taxes.values()):
                del taxes["Bacteria"]
            if "Archaea" in taxes and taxes["Archaea"] < .05 * sum(taxes.values()):
                del taxes["Archaea"]
            if "Bacteria" in taxes and taxes["Bacteria"] >= .95 * sum(taxes.values()):
                taxes = "Bacteria"
            elif "Archaea" in taxes and taxes["Archaea"] >= .95 * sum(taxes.values()):
                taxes = "Archaea"
            elif "Bacteria" in taxes and "Archaea" in taxes and taxes["Bacteria"] + taxes["Archaea"] >= .95 * sum(taxes.values()):
                taxes = f"Archaea {taxes['Archaea']}/ Bacteria {taxes['Bacteria']}"
            elif len(taxes) == 1:
                taxes = list(taxes.keys())[0]
                if taxes.startswith("cellular root_cellular organisms;"):
                    taxes = "; ".join(taxes.split("; ")[1:])
            else:
                mixed = "m"
                consensus = []
                for tax in taxes:
                    stax = tax.split("; ")
                    for i in range(0, len(stax)):
                        if "species_" in stax[i]:
                            stax[i] = stax[i].split("s_")[1]
                        try:
                            consensus[i]
                        except IndexError:
                            consensus.append({stax[i] : taxes[tax]})
                        else:
                            try:
                                consensus[i][stax[i]] += taxes[tax]
                            except:
                                consensus[i][stax[i]] = taxes[tax]
                        if not "_" in stax[i]:
                            break
                matched = []
                for level in consensus:
                    # level = set(level)
                    if len(level) == 1:
                        matched.append(list(level.keys())[0])
                    else:
                        hits_str = []
                        for hit in level:
                            hits_str.append(f"{hit}: {level[hit]}")
                        matched.append(", ".join(hits_str))
                        break
                if len(matched) > 2:
                    taxes = "; ".join(matched[1:])
                else:
                    taxes = "Poor resolution"
        seq2tax_dict[seq] = (perf, mixed, taxes)


    with open(seq2fulltax_file, "a") as out_fh:
        for seq, tax in seq2tax_dict.items():
            out_fh.write(f"{seq}\t{tax[0]}\t{tax[1]}\t{tax[2]}\n")
            seq2fulltax_dict[seq] = tax

