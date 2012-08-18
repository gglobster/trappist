# script to rename contigs in multifasta files

<<<<<<< HEAD
from genomes import all as genomes
from libs.common import load_multifasta, write_fasta

for genome in genomes:
    print genome['file']
    file_path = "data/genomes/"+genome['file']
    outfile_path = "data/renamed/"+genome['file']
    contigs = load_multifasta(file_path)
    renamed = []
    counter = 1
    for contig in contigs:
        contig.id = genome['name']+"_"+str(counter)
        contig_path = "data/contigs/"+contig.id+".fas"
        write_fasta(contig_path, contig)
        renamed.append(contig)
        counter +=1
    write_fasta(outfile_path, renamed)
=======
from sys import argv
from libs.common import load_multifasta, write_fasta, ensure_dir

from genomes import all as genome_list

origin_dir = "data/"+argv[1]+"/"
destin_dir = "data/"+argv[2]+"/"

ensure_dir([destin_dir])

for genome in genome_list:

    print genome['name'],

    while True:

        if genome['input'] == 'mfas':
            try:
                contigs = load_multifasta(origin_dir+genome['file'])
            except IOError:
                print "failed to load file"
                break

        else:
            print "input not appropriate"
            break

        renamed = []
        counter = 1
        for contig in contigs:
            contig.id = genome['name']+"_"+str(counter)
            renamed.append(contig)
            counter +=1

        write_fasta(destin_dir+genome['file'], renamed)

        print "OK"
        break
>>>>>>> massive update
