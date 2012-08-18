def batch_make_DBs():
    """Make Blast DBs from FastA files in a directory.

    This lists all files present in a directory, detects which are valid
    FastA files and creates Blast databases from their contents.

    """
    # get directory contents
    fas_dir = dirs['mfas_contigs_dir']
    db_dir = dirs['blast_db_dir']
    contents = listdir(fas_dir)
    g_names = []
    for item in contents:
        try:
            records = load_multifasta(item)
            assert len(records) > 0
        except IOError: print "\t"+item, "rejected (error opening file)"
        except AssertionError: print "\t"+item, "rejected (no FastA records)"
        except Exception: raise Exception
        else:
            print "\t"+item, "...",
            pattern = re.compile('^(.+)_.+')
            match = re.match(pattern, item)
            g_name = match.group(1)
            make_blastDB(db_dir+g_name, fas_dir+item, 'nucl')
            print "DB ready"
            g_names.append(g_name)
    return g_names