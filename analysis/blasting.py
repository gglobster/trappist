__author__ = 'GG'

def make_blastDB(db_path, genome_name, seq_file, db_type):
    """Create a BLASTable database if it doesn't already exist."""
    import sys, subprocess
    from other import ensure_dir
    dir_report = ensure_dir(db_path)
    if dir_report['status'] is 1:
        print dir_report['message']
        sys.exit()
    elif dir_report['status'] is 0: ### could use the absolute path! also, take out the hard-coding! 
        attempts = 0
        while attempts < 2:
            try: # check that the blast database files are all present        
                open('blast/'+genome_name+'.nhr')
                open('blast/'+genome_name+'.nin')
                open('blast/'+genome_name+'.nog')
                open('blast/'+genome_name+'.nsd')
                open('blast/'+genome_name+'.nsi')
                open('blast/'+genome_name+'.nsq')
            except IOError:
                attempts += 1
                status = 1
                cline = 'makeblastdb -in '+seq_file+' -dbtype '+db_type+' -title '\
                        +seq_file+' -out '+db_path \
                        +genome_name+' -parse_seqids'
                child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                         shell=True)
                output, error = child.communicate()
                message = {'output': output, 'error': error}
            else: 
                attempts += 2
                status = 0
                message = 'database exists'
    return {'message': message, 'status': status} 

def local_blastn(query_file, out_file, database, prefs):
    """Perform blastn on local database."""
    import subprocess
    from Bio.Blast.Applications import NcbiblastnCommandline
    cline = NcbiblastnCommandline(query=query_file, db=database,
                                  evalue=prefs['evalue'], out=out_file,
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    status = {'output': output, 'error': error}
    return status

def blast_record_set(db_name, fasta_records, blast_prefs):
    """Loop through fasta entries and blast against database."""
    from analysis.sequence_file_ops import write_fasta
    from analysis.blasting import local_blastn
    from analysis.text_manipulation import adaptive_list_load
    db_path = blast_prefs['db_path']
    evalue_pref = blast_prefs['evalue']
    score_pref = blast_prefs['score']
    len_pref = blast_prefs['length']
    matches = {}
    for query_record in fasta_records:
        # Blast each record against the current database.
        query_file = 'temp.fas'
        out_file = 'temp.blast'
        write_fasta(query_file,query_record)
        try:
            status = local_blastn(query_file, out_file, db_name, blast_prefs)
        except: raise
        else: # Parse BLAST results.
            results = adaptive_list_load(out_file,0,(2,3,4,9,10,11,12))
            for line in results:
                # Evaluate match details.
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
                    # Add match to dictionary with query id as key
                    if m_start > m_end: m_orient = '-'
                    else: m_orient = '+'
                    match = {'contig_id': contig_ID,
                              'match_p100': match_p100,
                             'length': length, 'm_orient': m_orient,
                             'evalue': evalue, 'bitscore': bitscore}
                    try: matches[query_record.id].append(match)
                    except KeyError: matches[query_record.id] = [match]
    return matches

