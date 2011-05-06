__author__ = 'GG'

def make_blastDB(db_path, genome_name, seq_file, db_type):
    """Create a BLASTable database if it doesn't already exist."""
    import os, subprocess
    from other import ensure_dir
    abs_path, dir_report = ensure_dir(db_path)
    dbfile_path = os.path.join(abs_path, genome_name)
    attempts = 0
    message, status = None, None
    while attempts < 2:
        try: # check that the blast database files are all present
            blast_db_ext_set = ['.nin', '.nog', '.nsd', '.nsi', '.nsq',
                                '.nhr']
            for blast_db_ext in blast_db_ext_set:
                open(dbfile_path+blast_db_ext)
        except IOError:
            attempts += 1
            status = 1
            cline = 'makeblastdb -in '+seq_file+' -dbtype '+db_type \
                    +' -title '+seq_file+' -out '+db_path \
                    +genome_name+' -parse_seqids'
            child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                     shell=True)
            output, error = child.communicate()
            message = {'output': output, 'error': error}
        else:
            attempts += 2
            status = 0
            message = 'database exists'
    report = {'message': message, 'status': status}
    return dbfile_path, report

def local_blastn(query_file, out_file, dbfile_path, prefs):
    """Perform blastn on local database."""
    import subprocess
    from Bio.Blast.Applications import NcbiblastnCommandline
    cline = NcbiblastnCommandline(query=query_file, db=dbfile_path,
                                  evalue=prefs['evalue'], out=out_file,
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    status = {'output': output, 'error': error}
    return status

def blast_record_set(dbfile_path, fasta_records, blast_prefs):
    """Loop through fasta entries and blast against database."""
    import os
    from analysis.sequence_file_ops import write_fasta
    from analysis.blasting import local_blastn
    matches = {}
    for query_record in fasta_records:
        query_file = 'temp.fas'
        out_file = 'temp.blast'
        write_fasta(query_file,query_record)
        try:
            status = local_blastn(query_file, out_file, dbfile_path,
                                  blast_prefs) 
        except: raise
        else:
            query_matches = parse_blast_out6(out_file,blast_prefs)
            matches[query_record.id] = query_matches
        finally:
            os.remove(query_file)
            os.remove(out_file)
    return matches

def parse_blast_out6(out_file, blast_prefs):
    """Parse BLAST results formatted in Tabular format (outfmt=6)."""
    from analysis.text_manipulation import adaptive_list_load
    score_pref = blast_prefs['score']
    len_pref = blast_prefs['length']
    query_matches = []
    results = adaptive_list_load(out_file, 0, (1, 2, 3, 8, 9, 10, 11))
    for line in results:
        contig_ID = line[0]
        match_p100 = float(line[1])
        length = int(line[2])
        m_start = int(line[3])
        m_end = int(line[4])
        evalue = line[5]
        bitscore = float(line[6])
        if bitscore < score_pref:
            pass # this is bad
        elif length < len_pref:
            pass # this is bad too
        else:
            if m_start > m_end: m_orient = '-'
            else: m_orient = '+'
            match = {'contig_id': contig_ID,
                      'match_p100': match_p100,
                     'length': length, 'm_orient': m_orient,
                     'evalue': evalue, 'bitscore': bitscore}
            query_matches.append(match)
    return query_matches