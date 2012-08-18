## script to combine several fasta sequences into a single one in a specific order

from sys import argv
from libs.common import load_fasta, write_genbank
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

origin_dir = "data/"+argv[1]+"/"
destin_file = origin_dir+argv[2]+".gbk"
base_name = argv[3]

# adapt this part
order = [(7,0), (17,1), (15,0), (16,1), (11,0), (6,0), (5,0), (1,0), (2,0), (3,0), (9,0), (13,0), (10,0), (4,0), (8,0), (12,0)]

filename = origin_dir+base_name+str(order[0][0])+".fas"
record = load_fasta(filename)

if order[0][1]:
    record = record.reverse_complement()
    c_note = '_RC'
else:
    c_note = ''

space_loc = FeatureLocation(0, len(record.seq))
quals = {'locus_tag': 'ctg_'+str(order[0][0])+c_note}
feature = SeqFeature(location=space_loc, type='contig',
                     id='ctg_'+str(order[0][0])+c_note, qualifiers=quals)
record.features.append(feature)

c_note = ''

for index in order[1:]:
    filename = origin_dir+base_name+str(index[0])+".fas"
    new_rec = load_fasta(filename)
    if index[1]:
        new_rec = new_rec.reverse_complement()
        c_note = '_RC'
    else:
        c_note = ''
    quals = {'locus_tag': 'ctg_'+str(index[0])+c_note}
    space_loc = FeatureLocation(len(record.seq),
                                len(record.seq)+len(new_rec.seq))
    feature = SeqFeature(location=space_loc, type='contig',
                         id='ctg_'+str(index[0])+c_note, qualifiers=quals)
    record.features.append(feature)

    record += new_rec
    record += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

record.id = argv[2]
record.name = argv[2]
record.seq.alphabet = generic_dna

write_genbank(destin_file, record)

