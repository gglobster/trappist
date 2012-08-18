import datetime
import sys,subprocess
from sys import argv
from Bio import SeqIO

def batch_iterator(iterator, batch_size) :
	# copied from http://biopython.org/wiki/Split_large_file
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def parcel_fasta(source_file, parcel_dir, nickname, parcel_size):
	# adapted from http://biopython.org/wiki/Split_large_file
	record_iter = SeqIO.parse(open(source_file),"fasta")
	parcel_files = []
	for i, batch in enumerate(batch_iterator(record_iter, parcel_size)) :
		filename = nickname+"_%i.fasta" % (i+1)
		handle = open(parcel_dir+filename, "w")
		count = SeqIO.write(batch, handle, "fasta")
		handle.close()
		parcel_files.append(filename)
		print "\twrote %i records to %s" % (count, filename)
	return parcel_files


# Variables

source_file = argv[1]		# must be multifasta
parcel_dir = argv[2]		
nickname = argv[3]

parcel_size = 10000

### Main script ###

print "--- Parceling to subsets ---"
print datetime.datetime.now()

parcel_files = parcel_fasta(source_file, parcel_dir, nickname, parcel_size)

print datetime.datetime.now()
print len(parcel_files), "parcel files generated"