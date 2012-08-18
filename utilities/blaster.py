# script to do batch blast (adapt blast flavor as required)

from os import path
from sys import argv
from libs.common import make_blastDB, local_tblastn_2file,\
    load_genbank, write_fasta, ensure_dir

from genomes import all as genome_list

data_dir = "data/"+argv[1]+"/"
genome_dir = "data/"+argv[2]+"/"
infile = data_dir+argv[3]

blast_out = data_dir+"blast_out/"

ensure_dir([blast_out])

# for record in list:
for genome in genome_list:
    print genome['name'],
    # mame a blast DB if it doesn't already exist
    genome_path = genome_dir+genome['file']
    dbfile_path = "data/blast_db/"+genome['name']
    if genome['input'] == 'cgbk':
        print "ignoring cgbk file"
    else:
        while True:
            if not path.exists(dbfile_path+".nhr"):
                if genome['input'] == 'gbk':
                    try:
                        print "converting,",
                        record = load_genbank(genome_path)
                    except IOError:
                        print "failed to load Genbank file"
                        break
                    else:
                        try:
                            genome_path = genome_dir+genome['name']+".fas"
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
                outfile = blast_out+genome['name']+".txt"
                prefs = {'evalue': 0.001, 'outfmt_pref': 6}
                print "blasting,",
                local_tblastn_2file(infile, dbfile_path, outfile, prefs)
            except Exception:
                print "failed to blast"
                break
            else:
                print "OK"
                break
