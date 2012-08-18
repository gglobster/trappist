#### Using EFetch to download full records from Entrez ####

from sys import argv
from Bio import Entrez, SeqIO

# Tell NCBI who we are
Entrez.email	= "Geraldine_VdAuwera@harvard.hms.edu"	

## Customizable bits ##
seqdir = argv[1]				# directory to save files to
infile = argv[2]				# input file with accession numbers
dbase = "nucleotide"			# database to search
rtype = "gb"					# record type

### Functions ###

## EnsureDir : Creates directory if it doesn't exist ##
# def EnsureDir(directory) :
#     dirhandle 	= os.path.dirname(directory)
#     if not os.path.exists(dirhandle):
#         os.makedirs(dirhandle)
        
## Strippa : Cleans up raw data from file read-in ##
def Strippa1(rawlines) :
	all_lines = []
	for line in rawlines:
		stripline = line.strip("\n")
		stripline = stripline.strip("\r")
		line_list = stripline.split("\t")
		all_lines.append(line_list[0])
	return all_lines

## EFetcher : Downloads sequence records by accession number ##
def EFetcher(rec_id,counter) :
	counter = counter + 1
	fname	= seqdir + rec_id + ".gbk"
	# fetch the record
	print "    " + str(counter) + ". downloading " + rec_id
	net_handle = Entrez.efetch(db=dbase, id=rec_id, rettype=rtype)
	out_handle = open(fname, "w")
	out_handle.write(net_handle.read())
	out_handle.close()
	net_handle.close()
	return counter
	print "        saved to file."

### Main program ###

print "### BatchFetcher ###"

# check that the save-to directory exists
#EnsureDir(seqdir) 

# confirm setup
print "Loading accession numbers from file " + infile + ", saving to " + seqdir

# load in the accession numbers
inlist		= open(infile, 'r')
rawlines	= inlist.readlines()
inlist.close()
rec_list	= Strippa1(rawlines)

# recap setlist
print "Found " + str(len(rec_list)) + " records to download :"
print rec_list

# initiate record count
counter		= 0

# run batch fetcher
for rec_id in rec_list :
	counter	= EFetcher(rec_id, counter) 

# confirm complete stop
print "BatchFetcher has downloaded " + str(counter) + " records to file."