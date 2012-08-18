# script to find the distance between query hit loci on a sequence

import re
from os import path
from sys import argv
from libs.common import make_blastDB, local_tblastn_2file, local_tblastx_2file, \
    local_blastn_2file, load_genbank, load_multifasta, write_fasta, \
    ensure_dir, from_dir, read_array, blast_dtypes

data_dir = "data/"+argv[1]+"/"
dir_in = data_dir+argv[2]+"/"
infile = data_dir+argv[3] # must be a fasta file with query sequences
file_ext = argv[4]
blast_mode = argv[5]

if len(argv) > 5:
    blast_mode = argv[5]
else:
    blast_mode = 'n' # nucleotide blast by default

blast_out = data_dir+"blast_out/"

ensure_dir([blast_out])

queries = load_multifasta(infile)

filenames = from_dir(dir_in, re.compile(r'.*\.'+file_ext))

for filename in filenames:

    rec_name = filename[:filename.find("."+file_ext)]
    print rec_name,

    genome_path = dir_in+filename
    dbfile_path = "data/blast_db/"+rec_name

    while True:
        if not path.exists(dbfile_path+".nhr"):
            if file_ext == 'gbk':
                try:
                    print "converting,",
                    record = load_genbank(genome_path)
                except IOError:
                    print "failed to load Genbank file"
                    break
                else:
                    try:
                        genome_path = dir_in+rec_name+".fas"
                        write_fasta(genome_path, record)
                    except Exception:
                        print "failed to write Fasta file"
                        break
            try:
                print "making a DB,",
                make_blastDB(dbfile_path, genome_path, 'nucl')
            except IOError:
                print "failed to make DB"
                break
        try:
            # blastx against each genome DB
            outfile = blast_out+rec_name+".txt"
            prefs = {'evalue': 0.001, 'outfmt_pref': 6}
            print "blasting,",
            if blast_mode == 'n':
                local_blastn_2file(infile, dbfile_path, outfile, prefs)
            elif blast_mode == 'tx':
                local_tblastx_2file(infile, dbfile_path, outfile, prefs)
            elif blast_mode == 'tn':
                local_tblastn_2file(infile, dbfile_path, outfile, prefs)
        except Exception:
            print "failed to blast"
            break
        else:
            print "finding range,"
            try:
                # load blast results from file
                rec_array = read_array(outfile, blast_dtypes)
            except IOError:
                print "failed to load blast results"
                break
            except Exception:
                print "failed to load blast results (unknown error)"
                break
            else:
                print rec_array[0]
                # take first line (= best hit) only
                try:
                    line = rec_array[0]
                except IndexError:
                    print "empty blast results"
                    break
                if line[8] < line[9]:
                    q_start, q_stop = line[8]-1, line[9]
                    rev_flag = False
                else:
                    q_start, q_stop = line[9]-1, line[8]
                    rev_flag = True
            break

# TODO: glomp results and output range