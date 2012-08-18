### Chop a sequence into fixed-sized segments
### Output as multifasta, then make a blast DB

import sys,subprocess
from sys import argv
from Bio import SeqIO, GenBank
from Bio.SeqRecord import SeqRecord

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
    cline = "makeblastdb -in "+ infile +" -dbtype "+ db_type +" -title "+ infile +" -out "+ name +" -parse_seqids"
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True) 
    output, error = child.communicate()
    return output 
    
# Get the input file name from the CLI
input_file = argv[1]

# Get the desired segment size from the CLI
seg_size = int(argv[2])

# Get the genome nickname from the CLI
genome_nick = argv[3]

print "chop sequence in ", argv[1], " to size ", argv[2]

# Create output file 
output_handle = input_file+'_chop_'+str(seg_size)+'.fas'
output_file = open(output_handle, 'w')

# Load the genome into memory
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
print len(seg_coords), "segments created"

# Extract the segment sequences, writing to file as we go
for coord_pair in seg_coords:
	start = coord_pair[0]
	stop = coord_pair[1]
	seg_seq = genome[0].seq[start:stop]
	seg_record = SeqRecord(seg_seq, id=genome_nick+"_"+str(start)+"_"+str(stop))
	#print seg_record.id
	SeqIO.write(seg_record, output_file, 'fasta')
print "written to", output_handle

# Make a BLAST DB from the output file
db_name = genome_nick+str(seg_size)
make_blastDB(db_name, output_handle, "nucl")
print db_name, "DB created"


