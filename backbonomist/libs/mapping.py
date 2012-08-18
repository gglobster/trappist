from os import listdir
from common import ensure_dir
from array_tetris import offset_q2r_coords
from drawing import contig_draw, pairwise_draw
from loaders import load_genbank
import numpy as np

def prep_maps(run_ref, run_id, timestamp, g_select, r_root_dir, run_dirs,
              genomes, fixed_dirs, segtype, min_size, fct_flags,
              fct_colors, idpt): 
    """Set up generation of various maps."""
    # set inputs and outputs
    ref_ctg_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    ref_gbk = run_ref.gbk
    cst_root = run_root+run_dirs['scaffolds_dir']+ref_ctg_n+"/"
    segments_root = run_root+run_dirs['aln_seg_dir']+ref_ctg_n+"/"
    ctg_segs_root = segments_root+"contigs/"
    cst_segs_root = segments_root+"constructs/"
    maps_root = run_root+run_dirs['maps_dir']+ref_ctg_n+"/"
    ctg_aln_maps_root = maps_root+"contig_alns/"
    cst_ann_maps_root = maps_root+"constructs_annot/"
    cst_aln_maps_root = maps_root+"constructs_aln/"
    ensure_dir([cst_root, ctg_segs_root, cst_segs_root, maps_root,
                ctg_aln_maps_root, cst_ann_maps_root, cst_aln_maps_root])
    print " ", ref_ctg_n, "...",
    # log
    logstring = "".join(["\n\n# Generate maps @", timestamp, "\n\n"])
    run_ref.log(logstring)
    # map of reference with segment details
    map_ref_segs(run_ref, run_id, r_root_dir, run_dirs, min_size,
                 fct_flags, fct_colors, idpt)
    # log
    logstring = "ref_map"
    run_ref.log(logstring)
    print logstring
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        # log
        logstring = "".join(["\n", g_name])
        run_ref.log(logstring)
        while True:
            try:
                if g_name in g_select: pass
                else: break
            except TypeError:
                pass
            print "\t", g_name, "...",
            scaff_gbk = cst_root+g_name+"_"+ref_ctg_n+"_scaffold.gbk"
            ctg_aln_maps_dir = ctg_aln_maps_root+g_name+"/"
            ensure_dir([ctg_aln_maps_dir])
            # maps of contigs aligned to reference
            logstring = "ctg_aln"
            print logstring,
            logstring = "".join(["\t", logstring])
            run_ref.log(logstring)
            map_ctg_alns(run_ref, ref_gbk, genome, ctg_segs_root,
                         ctg_aln_maps_dir, fixed_dirs, segtype, min_size,
                         fct_flags, fct_colors, idpt) 
            # map of scaffold construct
            logstring = "cst_ant"
            print logstring,
            logstring = "".join(["\t", logstring])
            run_ref.log(logstring)
            map_cst_annot(run_ref, genome, scaff_gbk, cst_ann_maps_root,
                          fct_flags, fct_colors)
            # map of construct aligned to reference
            logstring = "cst_aln"
            print logstring
            logstring = "".join(["\t", logstring])
            run_ref.log(logstring)
            map_cst_aln(run_ref, ref_gbk, genome, scaff_gbk, cst_segs_root,
                        cst_aln_maps_root, segtype, min_size, fct_flags,
                        fct_colors, idpt)
            break

def map_ctg_alns(run_ref, ref_gbk, genome, ctg_segs_root, maps_root,
                 fixed_dirs, segtype, min_size, fct_flags, fct_colors, idpt):
    """Generate maps of contigs aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = run_ref.name
    segs_root = ctg_segs_root+g_name+"/"
    ctgs_dir = fixed_dirs['gbk_contigs_dir']+g_name+"/"
    # list genbank files in matches directory
    try:
        dir_contents = listdir(segs_root)
    except OSError:
        msg = "\nWARNING: no matching segments"
        run_ref.log(msg)
        print msg
    else:
        for ctg_num in dir_contents:
            ctg_gbk = ctgs_dir+g_name+"_"+ctg_num+".gbk"
            seg_file = segs_root+ctg_num+"/"+ctg_num+"_"+ref_ctg_n+"_segs.txt"
            map_file = maps_root+g_name+"_"+ctg_num+"_vs_"+ref_ctg_n+".pdf"
            # start mapping
            try:
                # load segments TODO: add idp-based clumping
                segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
                # deactivate offsetting
                g_offset = (0,0)
                q_invert = False
                # generate graphical map
                pairwise_draw(ref_ctg_n, g_name+"_"+ctg_num, ref_gbk, ctg_gbk,
                             segdata, map_file, q_invert, g_offset, 'dual',
                             'dual', 'm', 'fct', 'fct', min_size,
                             fct_flags, fct_colors, idpt)
            except IOError:
                msg = "\nERROR: could not load segments data"
                run_ref.log(msg)
                print msg
            except StopIteration:
                msg = "\nERROR: could not make map"
                run_ref.log(msg)
                print msg

def map_cst_annot(run_ref, genome, scaff_gbk, maps_root, fct_flags,
                fct_colors):
    """Generate map of annotated scaffold construct."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = run_ref.name
    map_file = maps_root+g_name+"_<"+ref_ctg_n+".pdf"
    # start mapping
    try: open(scaff_gbk, 'r')
    except IOError:
        msg = "WARNING: No scaffold construct to map"
        run_ref.log(msg)
        print msg
    else:
        contig_draw(g_name, scaff_gbk, map_file, 'm', 'fct', fct_flags,
                fct_colors)

def map_cst_aln(run_ref, ref_gbk, genome, scaff_gbk, segs_root, maps_root,
                segtype, min_size, fct_flags, fct_colors, idpt):
    """Generate map of construct aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = run_ref.name
    seg_file = segs_root+g_name+"/"+g_name+"_"+ref_ctg_n+"_segs.txt"
    map_file = maps_root+g_name+"_vs_"+ref_ctg_n+".pdf"
    # start mapping
    try: open(scaff_gbk)
    except IOError:
        print "WARNING: No scaffold construct to map"
    else:
        try:
            # load segments TODO: add idp-based clumping
            segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
        except IOError:
                msg = "\nERROR: could not load segments data"
                run_ref.log(msg)
                print msg
        except StopIteration:
                msg = "\nERROR: could not make map"
                run_ref.log(msg)
                print msg
        else:
            # offset coordinates where desired
            try:
                g_offset = genome['offset']
                if g_offset[0] != 0 or g_offset[1] != 0:
                    q_len = len(load_genbank(scaff_gbk).seq)
                    segdata = offset_q2r_coords(segdata, q_len, g_offset,
                                                segtype)
                # determine whether to flip the query sequence (negative offset)
                if g_offset[1] < 0:
                    q_invert = True
                else:
                    q_invert = False
            except KeyError:
            	g_offset = (0,0)
                q_invert = False
            # generate graphical map
            pairwise_draw(ref_ctg_n, g_name, ref_gbk, scaff_gbk, segdata,
                         map_file, q_invert, g_offset, 'dual', 'dual', 'm',
                         'fct', 'fct', min_size, fct_flags, fct_colors, idpt)

def map_ref_segs(run_ref, run_id, r_root_dir, run_dirs, min_size,
                 fct_flags, fct_colors, idpt): 
    """Generate map of reference contig with segment details.

    This provides a comparison of the original reference and the
    re-annotated version.

    """
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    ori_file = run_ref.file
    ref_maps_root = run_root+run_dirs['ref_map_dir']
    ensure_dir([ref_maps_root])
    gbk_file = run_root+run_dirs['ref_gbk_dir']+ref_n+"_re-annot.gbk"
    map_file = ref_maps_root+ref_n+"_ref.pdf"
    # start mapping
    try:
        # make mock segment, full-length with 100% id
        record = load_genbank(gbk_file)
        length = len(record.seq)
        segdata = [[1, length, 1, length, 100]]
        # deactivate offsetting
        g_offset = (0,0)
        q_invert = False
        # generate graphical map
        pairwise_draw(ref_n+"_ra", ref_n+"_ori", gbk_file, ori_file,
                     segdata, map_file, q_invert, g_offset, 'dual', 'dual',
                     'm', 'fct', 'product', min_size, fct_flags,
                     fct_colors, idpt)
    except IOError:
        msg = "\nERROR: could not load segments data"
        run_ref.log(msg)
        print msg
    except StopIteration:
        msg = "\nERROR: could not make map"
        run_ref.log(msg)
        print msg
