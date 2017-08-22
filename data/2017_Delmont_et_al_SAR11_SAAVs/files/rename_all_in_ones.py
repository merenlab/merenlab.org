import os
import glob

all_in_ones = glob.glob("*.all_in_one")

for all_in_one in all_in_ones:
    gene_call = open(os.path.join(all_in_one, all_in_one.split(".")[0] + ".fasta"), "r").readline()
    gene_call = gene_call.strip().replace(">","")
    os.rename(all_in_one, gene_call + ".all_in_one")
