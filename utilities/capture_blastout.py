# script to capture sequences from the results of a batch blast

from sys import argv
from Bio.SeqRecord import SeqRecord
from libs.common import read_array, blast_dtypes, load_fasta, write_fasta, ensure_dir

data_dir = "data/"+argv[1]+"/"
main_in = data_dir+argv[2]+"_results.txt"
main_out = data_dir+argv[2]+"_ctxt.fas"
ctx_dir = data_dir+"context/"
capture_span = int(argv[3])

ensure_dir([ctx_dir])

records = []

rec_array = read_array(main_in, blast_dtypes)

descript = "Context of "+argv[2]+" ("+argv[3]+" bp either side)"

for line in rec_array:

    query = line[0]
    subject = line[1]

    print subject

    rev_flag = False
    if line[8] < line[9]:
        q_start, q_stop = line[8]-1, line[9]
        rev_flag = False
    else:
        q_start, q_stop = line[9]-1, line[8]
        rev_flag = True

    c_start, c_stop = q_start-capture_span, q_stop+capture_span

    master_seq = load_fasta("data/contigs_fas/"+subject+".fas")

    if c_start < 0:
        c_start = 0
    if c_stop > len(master_seq.seq):
        c_stop = len(master_seq.seq)

    seq_bit = master_seq[c_start:c_stop]

    if rev_flag:
        seq_bit = seq_bit.reverse_complement()
    record = SeqRecord(id=subject, seq=seq_bit.seq, description=descript)
    records.append(record)

    rec_file = ctx_dir+subject+"_"+query+"_ctxt.fas"
    write_fasta(rec_file, record)

write_fasta(main_out, records)

### TODO: modify bb tools to accept multifasta as sole input