__author__ = 'GG'

def seq_subset_load(infile, subset_mode, subset_args):
    """Load a subset of sequence segments from a sequence file."""
    from analysis.other import coord_chop, get_seq_subset_by_coords
    from analysis.seq_features import feat_collect, feature_coords
    from analysis.sequence_file_ops import load_multifasta, surefmt_load, \
        write_fasta
    from analysis.text_manipulation import adaptive_list_load
    if subset_mode is 'flatfile':
        # in this case the sequence file MUST be multifasta
        try: subset = load_multifasta(infile)
        except: raise
        else:
            print "set of", len(subset), "sequence segments"
            subset_file = infile
    else:
        # load the query single sequence file (convert format if necessary)
        try: seq_record = surefmt_load(infile, 'fasta', 'generic_dna')
        except: raise
        else: print "query sequence loaded from", infile
        # load or generate coordinate pairs for target segments
        if subset_mode is 'coordinates':
            try:
                coords_file, header, columns =  subset_args
                coords_list = adaptive_list_load(coords_file, header, columns)
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
        elif subset_mode is 'size':
            try:
                size = subset_args
                coords_list = coord_chop(len(seq_record.seq), size, subset_mode)
            except: raise
            else: print len(coords_list), "segments generated to fit", size
        else:
            print "ERROR: A mode MUST be specified."
        # collect subset of sequence segments using coords_list
        try: subset = get_seq_subset_by_coords(seq_record,coords_list)
        except: raise
        else: print "subset of", len(subset), "sequence segments"
        # save subset to multifasta file for later use or reference
        subset_file = seq_record.id+'_subset.fas'
        try: write_fasta(subset_file, subset)
        except: raise
        else: print "subset written to fasta file", subset_file
    return subset, subset_file

def genome_sets_load(input_file, input_prefs, db_path):
    """Load genome datasets listed in an input file."""
    import sys
    from classes.analysis_obj import GenomeSet
    from analysis.sequence_file_ops import ensure_fasta
    from analysis.text_manipulation import adaptive_list_load
    from analysis.blasting import make_blastDB
    header = input_prefs['header']
    columns = input_prefs['columns']
    genomes_list = adaptive_list_load(input_file, header, columns)
    print "prepping BLAST databases"
    genome_sets = []
    for line in genomes_list:
        genome_name = line[0]
        seq_file = line[1]
        try: db_infile = ensure_fasta(seq_file)
        except: raise
        else: print "genome FASTA sequence available in", db_infile
        DB_report = make_blastDB(db_path, genome_name, seq_file, 'nucl')
        if DB_report['status'] is 1:
            print genome_name, ":", DB_report['message']['error']
            sys.exit()
        elif DB_report['status'] is 0:
            print genome_name, ":", DB_report['message']
        new_genome_set = GenomeSet(db_infile, genome_name)
        genome_sets.append(new_genome_set)
    print "   ", len(genome_sets),"databases ready to search"
    return genome_sets