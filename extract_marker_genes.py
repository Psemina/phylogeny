#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fasta file conversion
from nucleic acid to amino acid
for phylogenetic analysis of LAGLIDADG type group II introns

created on May 15 22:25 2018
"""

import glob
import pickle
from Bio import SeqIO

'''
f = open("cog_id_list.txt", "r")
cog_id_list = []
for fline in f:
    fline = fline.rstrip()
    cog_id_list.append(fline)
f.close()
#print(cog_id_list)

s = open("./swissprot/uniprot_sprot.dat", "r")
#s = open("./test.txt", "r")
rf_cog_dict = {}
rf_ids = []
flag1, flag2 = False, False
for sline in s:
    sline = sline.rstrip()
    if sline[0:2]=="//":
        rf_ids = []
        cog_id = ""
    if sline[0:2]=="DR" and "COG" in sline:
        cog_id = sline.split(";")[1].replace(" ", "")
        if cog_id in cog_id_list:
            flag1 = True
    if sline[0:2]=="DR" and "RefSeq" in sline:
        rf_ids.append(sline.split(";")[1].replace(" ", ""))
        flag2 = True
    if flag1==True and flag2==True and "RefSeq" not in sline:
        for rf_id in rf_ids:
            rf_cog_dict.update({rf_id:cog_id})
        flag1, flag2 = False, False
s.close()

with open("rf_cog_dict.pickle", "wb") as g:
    pickle.dump(rf_cog_dict, g)

'''
with open("rf_cog_dict.pickle", "rb") as h:
    rf_cog_dict = pickle.load(h)

files = glob.glob("../../../gbff/rsBacteria_representative_180515_release87/*.gbff")
for afile in files:
    for seq_record in SeqIO.parse(afile, "genbank"):
        num_of_taxon = 0
        entry_id = ">{0}-{1}-".format(seq_record.id, afile.split("/")[5])
        for x in seq_record.features:
            if x.type=="source" and num_of_taxon==0:
                for y in x.qualifiers["db_xref"]:
                    if "taxon" in y:
                        entry_id += y.split(":")[1]
                        num_of_taxon += 1
            elif x.type=="source" and num_of_taxon!=0:
                for y in x.qualifiers["db_xref"]:
                    if "taxon" in y:
                        entry_id += "-{0}".format(y.split(":")[1])
                        num_of_taxon += 1
            elif x.type=="CDS":
                if "protein_id" in x.qualifiers.keys() and "translation" in x.qualifiers.keys():
                    for an_id in x.qualifiers["protein_id"]:
                        if an_id in rf_cog_dict.keys():
                            for prot_num, prot_sequence in enumerate(x.qualifiers["translation"]):
                                append_id = "-{0}-{1}-{2}-{3}".format(num_of_taxon, an_id, prot_num+1, rf_cog_dict[an_id])
                                id_print = "{0}-{1}".format(entry_id, append_id)
                            print("{0}\n{1}".format(id_print, prot_sequence))

