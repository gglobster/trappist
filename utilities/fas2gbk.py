## script to reverse orientation of a sequence, including all feature annotations

import re
from sys import argv
from backbonomist.libs.loaders import load_genbank
from backbonomist.libs.writers import write_genbank

infile = "data/"+argv[1]
outfile = "data/"+argv[2]

record = load_fasta(origin_dir+"/"+filename)

    # make a genbank file of the contig
    gbk_file = "".join([destin_dir, rec_name, ".gbk"])
    record.name = rec_name
    record.id = rec_name
    record.seq.alphabet = generic_dna
    write_genbank(gbk_file, record)

