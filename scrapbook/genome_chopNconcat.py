### Chop a set of genomes into fixed-sized segments
### Output as concatenated multifasta, then make a blast DB

import sys,subprocess
from sys import argv
import numpy
from Bio import SeqIO, GenBank
from Bio.SeqRecord import SeqRecord

def td_txt_file_load(filename,header):
    """Load raw info from tab-delimited text file."""
    infile = open(filename, 'r') 
    head_count = 0
    while head_count < header :
        infile.readline()
        head_count  = head_count +1
    rawlines = infile.readlines()
    infile.close()
    return rawlines

def singleton_list_load(filename,header):
    rawlines_list = td_txt_file_load(filename,header)
    singletons_list = []
    for rawline in rawlines_list :   
        line = strippa_clean(rawline)  
        if line[0] == "" :
            break
        singletons_list.append(line[0])
    return singletons_list

def strippa_clean(line):
    """Clean up raw data from file read-in, by line, into list array."""
    stripline = line.strip("\n")
    stripline = stripline.strip("\r")
    line_list = stripline.split("\t")
    return line_list

def load_agnostic(seqfile):
    """Load single-record file of unspecified format."""
    # present version is single-record fasta or genbank only
    while True:
        try: seq_record = load_fasta(seqfile)
        except: print "not single-record fasta"
        else: 
            rec_type = 'fasta'
            print "found a fasta file"
            return seq_record, rec_type
            break
        try: seq_record = load_genbank(seqfile)
        except: print "not genbank"
        else: 
            rec_type = 'genbank'
            print "found a genbank file"
            return seq_record, rec_type
            break
        raise Exception("Cannot open query file!")

def load_fasta(seqfile):
    """Load single-record Fasta file."""
    input_handle = open(seqfile, 'rU')
    fasta_record = SeqIO.read(input_handle, 'fasta')
    assert fasta_record.id 
    input_handle.close()
    return fasta_record

def load_genbank(seqfile):
    """Load single-record GenBank file."""
    parser = GenBank.FeatureParser()
    input_handle = open(seqfile, 'rU')
    gb_record = parser.parse(input_handle)
    assert gb_record.id 
    input_handle.close()
    return gb_record

def make_blastDB(name, infile, db_type):
    """Make BLAST database from FASTA input file.""" 
    cline = "makeblastdb -in "+ infile +" -dbtype "+ db_type +" -title "+ infile +" -out "+"data/blast_db/"+ name +" -parse_seqids"
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True) 
    output, error = child.communicate()
    return output 

    

# Variables

acc_ids = argv[1]
seg_size = int(argv[2])
working_dir = argv[3]
db_name = argv[4]

chop_dir = working_dir+'chopped/'

# Master concatenated file
concat_filename = chop_dir+'concat_genomes.fas'
concat_file = open(concat_filename, 'w')

# Load genome names
genomes_list = singleton_list_load(acc_ids, 0)

for g_name in genomes_list:

	print "chop genome", g_name, " to size ", argv[2]
	
	# Create output file 
	output_handle = chop_dir+g_name+'_chop_'+str(seg_size)+'.fas'
	output_file = open(output_handle, 'w')
	
	# Load the genome into memory
	input_file = working_dir+g_name+'.gbk'
	genome = load_agnostic(input_file)
	genome_length = len(genome[0].seq)
	
	# Create list of segment coordinates
	seg_coords = []
	index = 0
	while index < genome_length:
		incr_index = index+seg_size
		pair = (index, incr_index)
		seg_coords.append(pair)
		index = incr_index
	print "\t"+str(len(seg_coords)), "segments created"
	
	# Extract the segment sequences, writing to file (single & concat) as we go
	for coord_pair in seg_coords:
		start = coord_pair[0]
		stop = coord_pair[1]
		seg_seq = genome[0].seq[start:stop]
		seg_record = SeqRecord(seg_seq, id=g_name+"_"+str(start)+"_"+str(stop))
		#print seg_record.id
		SeqIO.write(seg_record, output_file, 'fasta')
		SeqIO.write(seg_record, concat_file, 'fasta')
	print "written to", output_handle
	output_file.close()

concat_file.close()
print "-- CONCAT FILE DONE --"

# Make a BLAST DB from the output file
db_name = db_name+str(seg_size)
make_blastDB(db_name, concat_filename, "nucl")
print "--", db_name, "DB created --"


