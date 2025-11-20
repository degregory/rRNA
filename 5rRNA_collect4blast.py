#!/bin/env python3

import os

seq2fulltax_file = '/mnt/c/Weekly_MiSeq/rRNA/seq2fulltax.tsv'

# get existing full tax
seq2fulltax_dict = {}
if os.path.isfile(seq2fulltax_file):
    with open(seq2fulltax_file, "r") as in_fh:
        for line in in_fh:
            line = line.strip().split("\t")
            seq2fulltax_dict[line[0]] = (line[1], line[2], line[3])
print(len(seq2fulltax_dict))

hit_dict = {}
samps_dict = {}

Species_dict_dict = {}
Per_Species_dict_dict = {}
all_species = {}
no_hits = {}

def check_chim(query, seqs):
    chimcheck = None
    front = 0
    back = 0
    length = len(query)
    for seq in seqs:
        if seq:
            if not query == seq[1]:
                ol_count = 0
                for i in range(0, min(length, len(seq[1]))):
                    if seq[1][i] == query[i]:
                        ol_count += 1
                    else:
                        break
                if ol_count > front:
                    front = ol_count
                ol_count = 0
                for i in range(1, min(length, len(seq[1]))):
                    if seq[1][-i] == query[-i]:
                        ol_count += 1
                    else:
                        break
                if ol_count > back:
                    back = ol_count
                if (front + back) > length:
                    return 1

    return chimcheck

u_count = 0

for file in os.listdir():
    if file.endswith("uc.fa"):
        samp_name = file.split(".")[0]
        Species_dict_dict[samp_name] = {}
        Per_Species_dict_dict[samp_name] = {}
        total = 0

        with open(file, "r") as in_fh, open(samp_name+".tsv", "w") as out_fh:
            seqs = in_fh.read().split(">")
            for i in range(0, len(seqs)):
                if seqs[i]:
                    seq = seqs[i].split("\n")
                    seqs[i] = (seq[0], ''.join(seq[1:]))
                    
            for seq in seqs:
                if not seq:
                    continue
                # seq = seq.split("\n")
                # seq[1] = ''.join(seq[1:])
                count = int(seq[0].split("size=")[1])
                total += count
                
                species = "No match"
                flag = ''
                mixed = ''
                try:
                    species = seq2fulltax_dict[seq[1]][2]
                    flag = seq2fulltax_dict[seq[1]][0]
                    mixed = seq2fulltax_dict[seq[1]][1]
                except:
                    pass
                
                if (not ('Archaea' in species or "Bacteria" in species)) and ";" in species and not flag:
                    if not check_chim(seq[1], seqs):
                        flag = "u"
                        u_count += 1
                
                
                
                out_fh.write(f"{flag}\t{mixed}\t{species}\t{seq[0]}\t{seq[1]}\n")
                if species == "No match":
                    no_hits[seq[1]] = 1
                if mixed:
                    ss = species.split("; ")
                    if ": " in ss[-1]:
                        species = "; ".join(ss[:-1])
                if '/ Bacteria' in species:
                    species = 'Archaea/Bacteria' 
                if flag and not ('Archaea' in species or "Bacteria" in species) and ";" in species:
                    species = f"{seq[1]}\t{species}"
                else:
                    species = f"\t{species}"
                try:
                    samps_dict[species]
                except:
                    samps_dict[species] = [[], []]
                try:
                    hit_dict[species]
                except:
                    hit_dict[species] = [0, 0]
                hit_dict[species][0] += count
                samps_dict[species][0].append(file)
                if flag == "p":
                    hit_dict[species][1] += count
                    samps_dict[species][1].append(file)
                    try:
                        Per_Species_dict_dict[samp_name][species] += count
                    except:
                        Per_Species_dict_dict[samp_name][species] = count
                all_species[species] = 1
                try:
                    Species_dict_dict[samp_name][species] += count
                except:
                    Species_dict_dict[samp_name][species] = count
        Species_dict_dict[samp_name]['total'] = total

print(u_count)

if all_species:
    with open('rRNA_Species.tsv', "w") as rs_fh:
        rs_fh.write("Sequence\tTaxonomy\tTax Level\tPerfect hits\tpos samps\tper pos samps")
        rs_fh.write("\tpos reads\tper pos reads")

        for samps in Species_dict_dict:
            rs_fh.write(f"\t{samps}({Species_dict_dict[samps]['total']})\t")
        rs_fh.write("\n")
        rs_fh.write("\t\t\t\t\t\t\t")
        for samps in Species_dict_dict:
            rs_fh.write(f"\tCount\tAbundance")
        rs_fh.write("\n")
        sorted_species = sorted(all_species, key=lambda x: " ".join(x.split("\t")[::-1]))
        for spec in sorted_species:
            level = spec.split("\t")[-1].split('; ')
            if len(level) > 1 and not "_" in level[-1]:
                level = "Species"
            elif "species_" in spec:
                level = "Species"
            else:
                level = level[-1].split("_")[0]
            if level == 'clade':
                cc = 0
                for entry in spec.split('; ')[::-1]:
                    if "clade_" in entry:
                        cc += 1
                    else:
                        level = entry.split("_")[0].strip()+f"({cc})"
                        break
            rs_fh.write(f"{spec}\t{level}")
            
            perfect = 0
            if hit_dict[spec][1] + hit_dict[spec][1] > 0:
                perfect = 1
            rs_fh.write(f"\t{perfect}\t{len(set(samps_dict[spec][0]))}\t{len(set(samps_dict[spec][1]))}")
            rs_fh.write(f"\t{hit_dict[spec][0]}\t{hit_dict[spec][1]}")
            for samps in Species_dict_dict:
                per_hits = 0
                try:
                    per_hits = Per_Species_dict_dict[samps][spec]
                except:
                    pass
                try:
                    rs_fh.write(f"\t{Species_dict_dict[samps][spec]} ({per_hits})\t{Species_dict_dict[samps][spec]/Species_dict_dict[samps]['total']}")
                except:
                    rs_fh.write(f"\t\t")
            rs_fh.write("\n")
