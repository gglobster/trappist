#### Using EFetch to download full records from Entrez ####

# TODO: clean up & standardize to NGSU coding practices

from sys import argv
from Bio import Entrez
from libs.common import ensure_dir, load_genbank, write_fasta

# Tell NCBI who we are
Entrez.email	= "Geraldine_VdAuwera@harvard.hms.edu"

## Customizable bits ##
data_dir = "data/"+argv[1]+"/"    # directory to save files to
infile = data_dir+argv[2]			# input file with accession numbers
dbase = "nucleotide"			# database to search
rtype = "gb"					# record type

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
def EFetcher(rec_id, seqdir) :
    fname = seqdir+rec_id+".gbk"
    # fetch the record
    print ".",
    net_handle = Entrez.efetch(db=dbase, id=rec_id, rettype=rtype,
                               retmode = 'text') # retmode defaults to XML
    out_handle = open(fname, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    return fname

### Main program ###

print "### BatchFetcher ###"

# check that the save-to directory exists
#EnsureDir(seqdir)

# confirm setup
print "Loading accession numbers from file", infile, ", saving to ", data_dir

# load in the accession numbers
inlist = open(infile, 'r')
rawlines = inlist.readlines()
inlist.close()
rec_list = Strippa1(rawlines)

# recap setlist
print "Found " + str(len(rec_list)) + " records to download"

# initiate record count
counter = 0

# run batch fetcher
for rec_id in rec_list :

    print rec_id,
    counter += 1

    while True:
        try:
            fname = EFetcher(rec_id, data_dir)
        except Exception:
            print "Error retrieving record"
            break
        else:
            if rec_id[0:2] == 'NZ': # disposition for WGS record sets

                print "fetching WGS dataset",

                # create a dedicated directory
                seqdir = data_dir+rec_id+"/"
                ensure_dir([seqdir])

                # open genome record stub to get the contig count
                fname = data_dir+rec_id+".gbk"
                try:
                    stub = load_genbank(fname)
                except IOError:
                    print "Error loading", fname
                    break

                base_code = stub.annotations['wgs'][0][:10] # 7 if not NZ_
                ctg_num = int(stub.annotations['wgs'][-1][10:]) # 7

                records = []

                # fetch contig records
                ctg_count = 0
                while ctg_count < ctg_num:
                    # TODO: better formatting
                    ctg_count += 1
                    if ctg_count < 10:
                        ctg_id = base_code+'0000'+str(ctg_count)
                    elif ctg_count < 100:
                        ctg_id = base_code+'000'+str(ctg_count)
                    elif ctg_count < 1000:
                        ctg_id = base_code+'00'+str(ctg_count)
                    else: # shouldn't happen but hey...
                        ctg_id = base_code+'0'+str(ctg_count)
                    # fetch contig record
                    try:
                        fname = EFetcher(ctg_id[3:], seqdir) # 3 if not NZ_
                    except Exception:
                        print "Error retrieving record"
                    else:
                        try:
                            records.append(load_genbank(fname))
                        except Exception:
                            print "Error loading record"

                write_fasta(data_dir+rec_id+".fas", records)

            print "OK"
            break

# confirm complete stop
print "BatchFetcher has downloaded " + str(counter) + " records to file."