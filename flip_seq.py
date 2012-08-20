## script to reverse orientation of a sequence, including all feature annotations

from sys import argv
from backbonomist.libs.loaders import load_genbank
from backbonomist.libs.writers import write_genbank

infile = "data/"+argv[1]
outfile = "data/"+argv[2]

record = load_genbank(infile)

# make a genbank file of the contig
new_record = record.reverse_complement()
write_genbank(outfile, new_record)

