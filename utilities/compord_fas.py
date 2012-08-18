## script to combine several fasta sequences into a single one in a specific order

from sys import argv
from libs.common import  load_fasta, write_fasta

origin_dir = "data/"+argv[1]+"/"
destin_file = origin_dir+argv[2]+".fas"
base_name = argv[3]

# adapt this part
order = [(22, 0), (4, 0), (57, 1), (43, 1), (64, 0), (18, 0), (54, 0), (36, 1), (20, 1), (2, 1), (40, 1), (17, 1), (35, 1), (38, 1), (37, 1), (55, 1), (19, 1), (47, 1), (11, 0), (46, 0), (61, 0), (41, 1), (15, 0), (1, 1), (5, 1), (6, 0), (13, 1), (8, 0), (23, 0), (16, 1), (10, 0), (60, 0), (14, 0), (42, 0), (39, 0), (48, 0), (9, 1), (21, 0), (3, 1), (58, 1), (32, 0)]

filename = origin_dir+base_name+str(order[0][0])+".fas"
record = load_fasta(filename)
if order[0][1]:
    record = record.reverse_complement()

for index in order[1:]:
    filename = origin_dir+base_name+str(index[0])+".fas"
    new_rec = load_fasta(filename)
    if index[1]:
        new_rec = new_rec.reverse_complement()
    record += new_rec
    record += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

record.id = argv[2]

write_fasta(destin_file, record)