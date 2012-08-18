# script to translate sequences in multifasta files into proteins

from sys import argv
from libs.common import load_multifasta, write_fasta
from Bio.SeqRecord import SeqRecord

origin_dir = "data/"+argv[1]+"/"
in_file = origin_dir+argv[2]
outfile = in_file[:-4]+"_aa.fas"

proteins = []

for record in load_multifasta(in_file):
    aa_rec = SeqRecord(id=record.id, seq=record.seq.translate())
    proteins.append(aa_rec)

write_fasta(outfile, proteins)
