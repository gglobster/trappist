## script to combine several fasta files into a single one

import re
from sys import argv
from libs.common import from_dir, load_fasta, load_genbank, write_fasta

origin_dir = "data/"+argv[1]
destin_file = origin_dir+"/"+argv[2]+".fas"
file_ext = argv[3]

filenames = from_dir(origin_dir, re.compile(r'.*\.'+file_ext))

records = []

for filename in filenames:
    # load record
    if file_ext == 'fas':
        records.append(load_fasta(origin_dir+"/"+filename))
    elif file_ext == 'gbk':
        records.append(load_genbank(origin_dir+"/"+filename))

    print filename

write_fasta(destin_file, records)