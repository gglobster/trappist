### Blast Illumina data against a reference DB
### Filter out positive matches to separate file

import datetime
import sys,subprocess
from sys import argv
from Bio import SeqIO, GenBank
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Blast.Applications import NcbiblastnCommandline


def local_blastn_live(query_file, dbfile_path, prefs):
    """Perform tblastn on local database."""
    cline = NcbiblastnCommandline(query=query_file, db=dbfile_path,
                                  evalue=prefs['evalue'], outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def write_temp_query(temp_infile, stuff):
	# write out seq to temp text file 
	# not using Fasta to avoid extra overhead
	query = open(temp_infile, "w")
	query.write(stuff)
	query.close()

# Variables

source_file = argv[1]		# must be FASTQ
nickname = argv[2]			# avoid numbers in nickname
db_name = argv[3]

temp_infile = "data/temp/seq_in.txt"
dbfile_path = "data/blast_db/"+db_name
blast_out_path = "data/blast_out/"
match_file = blast_out_path+nickname+"_"+db_name+"_positive.txt"
no_match_file = blast_out_path+nickname+"_"+db_name+"_negative.txt"

blast_prefs = {'evalue': 0.001, 'outfmt_pref': 6, 'score': 100, 'length': 60}

# Blast each read against reference genome DB

read_count = 0
match_count = 0
no_match_count = 0
matches = open(match_file, "w")
no_match = open(no_match_file, "w")

print "--- Starting analysis run ---"
print datetime.datetime.now()

for title, seq, qual in FastqGeneralIterator(open(source_file)) :
	id_string = nickname+"_"+str(read_count)
	write_temp_query(temp_infile, ">"+id_string+"\n"+seq)
	output = local_blastn_live(temp_infile, dbfile_path, blast_prefs)
	if output == "":
		no_match.write(id_string+"\n")
		no_match_count +=1
	else:
		matches.write(output)
		match_count +=1
	read_count +=1
	if read_count%10==0:
		print read_count, "reads processed"
	if read_count == 500:
		break
matches.close()
no_match.close()

print "--- Finished analysis run ---"
print datetime.datetime.now()
print read_count, "reads blasted"
print match_count, "reads matched positive"
print no_match_count, "reads didn't match"

    