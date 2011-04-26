from shared.blasting import make_blastDB

__author__ = 'GG'

def seq_subset_load(infile, subset_mode, subset_args):
    """Load a subset of sequence segments from a sequence file."""
    import sys
    from shared.other import coord_chop, get_seq_subset_by_coords
    from shared.seq_features import feat_collect, feature_coords
    from shared.sequence_file_ops import surefmt_load, write_fasta
    from shared.text_manipulation import adaptive_list_load
    # load the query sequence file (convert format if necessary)
    try: seq_record = surefmt_load(infile, 'fasta', 'generic_dna')
    except: raise   
    else: print "query sequence loaded from", infile
    # load or generate coordinate pairs for target segments
    if subset_mode is 'file':
        try: 
            infile, header, columns =  subset_args
            coords_list = adaptive_list_load(infile, header, columns)
        except: raise
        else: print len(coords_list), "segments loaded from", infile
    elif subset_mode is 'features':
        try:
            feat_mode = subset_args
            features = feat_collect(infile, feat_mode)
            coords_list = feature_coords(features)
            print coords_list
        except: raise
        else: print len(coords_list),"features loaded from", infile
    else:
        try: 
            size = subset_args
            coords_list = coord_chop(len(seq_record.seq), size, subset_mode)
        except: raise
        else: print len(coords_list), "segments generated to fit", size
    # collect subset of sequence segments using coords_list
    try: subset = get_seq_subset_by_coords(seq_record,coords_list)
    except: raise
    else: print "subset of", len(subset), "sequence segments"
    # save subset to multifasta file for later use or reference
    subset_file = seq_record.id+'_subset.fas'
    try: write_fasta(subset_file, subset)
    except: raise
    else: print "subset written to fasta file", subset_file
    # clean up and move on
    return subset, subset_file

def genome_set_load(input_file, input_prefs, db_type):
    """Load a set of genomes listed in an input file."""
    import sys
    from classes.analysis_obj import GenomeSet
    from shared.sequence_file_ops import ensure_fasta
    from shared.text_manipulation import adaptive_list_load
    from shared.blasting import make_blastDB
    header = input_prefs['header']
    columns = input_prefs['columns']
    genomes_list = adaptive_list_load(input_file, header, columns)
    print "prepping BLAST databases"
    genome_sets = []
    for line in genomes_list:
        genome_name = line[0]
        seq_file = line[1]
        db_infile = ensure_fasta(seq_file)
        DB_report = make_blastDB(genome_name, seq_file, db_type)
        if DB_report['status'] is 1:
            print genome_name, ":", DB_report['message']['error']
            sys.exit()
        elif DB_report['status'] is 0:
            print genome_name, ":", DB_report['message']
        """Instantiate a new RamenCup object from line content."""
        new_ramen_cup   = RamenCup(data_dir,db_infile,cup_name)
        ramen_cups_RA.append(new_ramen_cup)
    print "   ", len(ramen_cups_RA),"databases ready to search"
    return ramen_cups_RA