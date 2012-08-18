### Blast batches of files with common root name ###

import os
import datetime
import sys,subprocess
from sys import argv
from Bio.Blast.Applications import NcbiblastnCommandline


def local_blastn_2file(query_file, dbfile_path, outfile, prefs):
    """Perform tblastn on local database."""
    cline = NcbiblastnCommandline(query=query_file, db=dbfile_path, out=outfile,
                                  evalue=prefs['evalue'], outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait for blastn to finish
    

# Variables

parcel_dir = argv[1]		# location of the parcel files
species_dir = argv[2]
nickname = argv[3]
db_name = argv[4]

dbfile_path = "data/blast_db/"+db_name
blast_out_path = "data/blast_out/"+species_dir+nickname+"/"

blast_prefs = {'evalue': 0.001, 'outfmt_pref': 6, 'score': 100, 'length': 60}

### Main script ### 

print "--- Starting batch BLAST run ---"
print datetime.datetime.now()

index = 1
while os.path.isfile(parcel_dir+nickname+"_"+str(index)+".fasta"):
	query_file = parcel_dir+nickname+"_"+str(index)+".fasta"
	outfile = blast_out_path+nickname+"_"+str(index)+"_blast.out"
	print "\tblasting", query_file
	local_blastn_2file(query_file, dbfile_path, outfile, blast_prefs)
	index +=1

print "--- Finished BLAST run ---"
print datetime.datetime.now()
print index, "parcel files blasted"

    