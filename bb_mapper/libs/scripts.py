from classes import Noodle
from aligning import mauve_pw_align
from array_tetris import process_segdata
from drawing import pairwise_draw, multi_draw

def align_pairwise(genomes, new_align, r_root_dir, g_root_dir,dirs, run, max_size, chop_mode, mauve_exec, mtype, segtype, min_size):
    """Make a pairwise alignment."""
    # set up directories
    aln_dir = r_root_dir+run+dirs['aln_segs']
    map_dir = r_root_dir+run+dirs['maps']
    # load inputs and process genomes
    seq_dir = g_root_dir
    ref = [Noodle(genome, seq_dir) for genome in genomes
           if genome['order'] == 1][0]
    query = [Noodle(genome, seq_dir) for genome in genomes
             if genome['order'] == 2][0]
    # align if needed
    print "Aligning", ref.name, "and", query.name, "...",
    if new_align:
        mauve_pw_align(ref, query, r_root_dir, g_root_dir, dirs, run, max_size, chop_mode, mauve_exec, mtype)
    # process data segments
    print "Processing segments ..",
    seg_file = aln_dir+ref.name+"_"+query.name+"_segs.txt"
    pair_data = process_segdata(seg_file, ref, query, segtype)
    print "OK"
    # map of query aligned to reference
    print "Mapping ...",
    map_file = map_dir+run+"_"+ref.name+"_vs_"+query.name+".pdf"
    pairwise_draw(ref, query, pair_data, map_file, 'dual', 'dual', 'm',
                  'fct', 'fct', idpt, fct_flags, fct_colors, min_size)
    print "OK\n"

def align_multi(genomes, new_align, r_root_dir, g_root_dir, dirs, run, max_size, chop_mode, mauve_exec, mtype, segtype, idpt, fct_flags, fct_colors, min_size):
    """Make a multiple alignment."""
    # set up directories
    aln_dir = r_root_dir+run+dirs['aln_segs']
    map_dir = r_root_dir+run+dirs['maps']
    # load inputs, process and pair up genomes
    seq_dir = g_root_dir
    g_pairs = []
    counter = 1
    while counter < len(genomes):
        g_pairs.append(([Noodle(genome, seq_dir) for genome in genomes if
                            genome['order'] == counter][0],
                        [Noodle(genome, seq_dir) for genome in genomes if
                            genome['order'] == counter+1][0]))
        counter +=1
    # process pairs
    for (ref, query) in g_pairs:
        # align if needed
        print "Aligning", ref.name, "and", query.name, "...",
        if new_align:
            mauve_pw_align(ref, query, r_root_dir, g_root_dir, dirs, run, max_size, chop_mode,
                           mauve_exec, mtype)
    # traverse genome pairs
    segdata_list = []
    counter = 0
    for (ref, query) in g_pairs:
        counter +=1
        # process data segments
        print "Processing pair", counter, "segments ...",
        p_seg_file = aln_dir+ref.name+"_"+query.name+"_segs.txt"
        pair_data = process_segdata(p_seg_file, ref, query, segtype)
        segdata_list.append(pair_data)
        print "OK"
        # make a pairwise map while we're at it
        print "Mapping pair", counter, "...",
        p_map_file = map_dir+run+"_"+ref.name+"_vs_"+query.name+".pdf"
        pairwise_draw(ref, query, pair_data, p_map_file, 'dual', 'dual', 'm',
                      'fct', 'fct', idpt, fct_flags, fct_colors, min_size)
        print "OK"
    # map of query aligned to reference
    print "Mapping multiple alignment...",
    map_file = map_dir+run+".pdf"
    multi_draw(g_pairs, segdata_list, map_file, idpt, fct_flags, fct_colors, min_size)
    print "OK\n"














































