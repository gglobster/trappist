## script to unpack contigs

from sys import argv
from libs.common import ensure_dir, load_multifasta, load_fasta,load_genbank, \
    write_fasta

from genomes import all as genome_list

origin_dir = "data/"+argv[1]+"/"
destin_dir = "data/"+argv[2]+"/"

ensure_dir([destin_dir])

for genome in genome_list:

    print genome['name'],

    origin_file = origin_dir+genome['file']

    while True:

        if genome['input'] == 'cgbk':
            print "ignoring cgbk file"
            break

        elif genome['input'] == 'gbk':
            try:
                records = [load_genbank(origin_file)]
            except IOError:
                print "failed to load file"
                break

        elif genome['input'] == 'fas':
            try:
                records = [load_fasta(origin_file)]
            except IOError:
                print "failed to load file"
                break

        elif genome['input'] == 'mfas':
            try:
                records = load_multifasta(origin_file)
            except IOError:
                print "failed to load file"
                break

        else:
            print "input not recognized"
            break

        for record in records:
            try:
                write_fasta(destin_dir+record.id+".fas", record)
            except Exception:
                print "failed to write contig file"
                break
            else:
                print record.id,

        print "OK"
        break


