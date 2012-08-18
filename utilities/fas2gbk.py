## script to convert sets of fasta files to genbank format

import re
from sys import argv
from libs.common import from_dir, ensure_dir, load_fasta, write_genbank
from Bio.Alphabet import generic_dna

origin_dir = "data/"+argv[1]
destin_dir = "data/"+argv[2]+"/"

ensure_dir([destin_dir])

filenames = from_dir(origin_dir, re.compile(r'.*\.fas.*'))

for filename in filenames:
    rec_name = filename[:filename.find('.fas')]
    record = load_fasta(origin_dir+"/"+filename)

    # make a genbank file of the contig
    gbk_file = "".join([destin_dir, rec_name, ".gbk"])
    record.name = rec_name
    record.id = rec_name
    record.seq.alphabet = generic_dna
    write_genbank(gbk_file, record)

