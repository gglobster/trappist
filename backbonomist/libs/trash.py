def unpack_gbk_concat_contigs(genome):
    """Unpack contigs concatenated with a sequence separator.

    This unpacking method consists of two main steps. The first is to
    detect occurrences of the sequence separator and collect the start and
    stop coordinates of each contig. The second is to use each pair of
    coordinates to extract the contig sequence and create a SeqRecord for
    that contig. Each contig record is then written to the master multifasta
    file and to its own single-record GenBank file.

    """
    # set up inputs
    infile = genome['file'] #TODO: make automated input loader (upstream)
    inpath = dirs['ori_g_dir']+infile
    # load in genome data
    genome_rec = load_genbank(inpath)
    g_string = genome_rec.seq
    # find split coordinates
    coord_pairs = multisplit_finder(g_string, separator)
    # prep output destinations
    mfas_dir = dirs['mfas_contigs_dir']
    gbk_dir = dirs['gbk_contigs_dir']+genome_rec.id+"/"
    ensure_dir(mfas_dir)
    ensure_dir(gbk_dir)
    fas_file = mfas_dir+genome_rec.id+"_contigs.fas"
    # split record
    counter = 0
    records = []
    for (start, stop) in coord_pairs:
        counter +=1
        new_record = genome_rec[start:stop]
        new_record.id = genome_rec.id+"_"+str(counter)
        records.append(new_record)  # for multifasta output, to blastdb
        gbk_file = gbk_dir+genome_rec.id+"_"+str(counter)\
                                        +"_"+str(start)\
                                        +"_"+str(stop)\
                                        +".gbk"
        write_genbank(gbk_file, new_record)
    write_fasta(fas_file, records)
    return genome_rec.id