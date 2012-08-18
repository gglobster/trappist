# script to process the results of a batch blast

import re
from os import path
from sys import argv
from Bio.SeqRecord import SeqRecord
from libs.common import from_dir, read_array, blast_dtypes, load_fasta, write_fasta

data_dir = "data/"+argv[1]+"/"
blast_out_dir = "data/"+argv[1]+"/blast_out/"
idp = int(argv[2])

main_out = open(data_dir+"comp_results.txt", 'w')

records_dict = {}

# list files in blast results directory
filenames = from_dir(blast_out_dir, re.compile(r'.*\.txt.*'))
for filename in filenames:
    counter = 0
    # load text
    rec_array = read_array(blast_out_dir+filename, blast_dtypes)
    # parse lines
    for line in rec_array:
        # if idp is higher than spec'd:
        if line[2] > idp:
            query = line[0]
            subject = line[1]
            # write line to compiled results file
            main_out.write("\t".join([str(item) for item in line])+"\n")
            outfile = data_dir+query+"_results.txt"
            if not path.exists(outfile):
                # create file
                out_handle = open(outfile, 'w')
            else:
                counter +=1
                out_handle = open(outfile, 'a')
            out_handle.write("\t".join([str(item) for item in line])+"\n")
            # extract sequence to array
            rev_flag = False
            if line[8] < line[9]:
                q_start, q_stop = line[8]-1, line[9]
                rev_flag = False
            else:
                q_start, q_stop = line[9]-1, line[8]
                rev_flag = True
            master_seq = load_fasta("data/contigs_fas/"+subject+".fas")
            seq_bit = master_seq[q_start:q_stop]
            if rev_flag:
                seq_bit = seq_bit.reverse_complement()
            record = SeqRecord(id=subject+"_"+str(counter), seq=seq_bit.seq)
            if query not in records_dict.keys():
                records_dict[query] = [record]
            else:
                records_dict[query].append(record)
# write out sequences
for query in records_dict.keys():
    seqfile_nt = data_dir+query+"_nt.fas"
    write_fasta(seqfile_nt, records_dict[query])



