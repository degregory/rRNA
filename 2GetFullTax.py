#!/usr/bin/env python
# coding: utf-8

import os
import sys
from Bio import Entrez
import time
Entrez.email = 'gregoryde@umsystem.edu'

taxid2tax_file = '/mnt/c/Weekly_MiSeq/rRNA/taxid2fulltax.tsv'


tax_ids = []
with open("blastn.tsv", "r") as in_fh:
    for line in in_fh:
        line = line.split("\t")
        tax_ids.append(line[9])

tax_ids = set(tax_ids)
print(len(tax_ids))
# get local taxid and associated full tax
fulltax_dict = {}
if os.path.isfile(taxid2tax_file):
    with open(taxid2tax_file, "r") as in_fh:
        for line in in_fh:
            line = line.split("\t")
            fulltax_dict[line[0]] = line[1]

# get full tax for new ids
tax_ids = [tax_id for tax_id in tax_ids if not tax_id in fulltax_dict]

print(len(tax_ids))

if tax_ids:
    with open(taxid2tax_file, "a") as out_fh:
        batch_size = 2500
        batches = [tax_ids[i: i + batch_size] for i in range(0, len(tax_ids), batch_size)]
        print(len(batches))
        count = 0
        
        for batch in batches:
            handle = Entrez.efetch(db='taxonomy', id=' '.join(batch), retmode='XML') # , api_key=api_key)
            ncbi_results = Entrez.read(handle)
            handle.close()
            time.sleep(7)
            for result in ncbi_results:
                full_tax = [(level['Rank'], level['ScientificName']) for level in result['LineageEx'] if level['Rank']]
                full_tax = "; ".join(["_".join(entry) for entry in full_tax]) + f"; {result['ScientificName']}"
                full_tax = full_tax.replace("no rank_", "clade_")
                # fulltax_dict[result['TaxId']] = full_tax
                out_fh.write(f"{result['TaxId']}\t{full_tax}\n")
