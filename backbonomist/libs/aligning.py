import os, subprocess, re
from os import listdir
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from common import ensure_dir
from parsing import mauver_load2_k0, parse_clustal_idstars
from array_tetris import chop_rows
from loaders import load_genbank, load_fasta
from writers import write_fasta

def align_clustal(file_name):
    """Make external call to ClustalW aligner."""
    cline = ClustalwCommandline(infile=file_name)
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse ClustalW errors
    return report

def align_muscle(infile_name, outfile_name, log_file):
    """Make external call to Muscle aligner."""
    cline = MuscleCommandline(input=infile_name, out=outfile_name, clw=True,
                              loga=log_file, quiet='y')
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse MUSCLE errors
    return report

def align_mauve(file_list, output, mauve_exec):
    """Make external call to Mauve aligner."""
    input_files = ' '.join(file_list)
    cline = mauve_exec+" --output="+ output +" "+input_files
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse Mauve errors
    return report

def align_ctg2ref(run_ref, run_id, timestamp, r_root_dir, run_dirs,
                  genomes, mauve_exec, max_size, chop_mode, mtype):
    """Align contigs pairwise to the reference contig."""
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    ref_ctg_file = run_ref.file
    mauve_root = run_root+run_dirs['mauve_out_dir']+ref_n+"/contigs/"
    segments_root = run_root+run_dirs['aln_seg_dir']+ref_n+"/contigs/"
    q_ctgs_root = run_root+run_dirs['match_out_dir']+ref_n+"/"
    ensure_dir([segments_root])
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Align contigs to ref @", timestamp, "\n"])
    run_ref.log(logstring)
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        ctgs_fas_dir = q_ctgs_root+g_name+"/"
        mauve_dir = mauve_root+g_name+"/"
        aln_segs_root = segments_root+g_name+"/"
        ensure_dir([mauve_dir])
        print "\t", g_name, "...",
        # log
        logstring = "".join(["\n", g_name])
        run_ref.log(logstring)
        # list genbank files in matches directory
        dir_contents = listdir(ctgs_fas_dir)
        for item in dir_contents:
            pattern = re.compile(r'.*_(\d*)\.fas$')
            match = pattern.match(item)
            if match:
                ctg_num = match.group(1)
                print ctg_num,
                logstring = "".join(["\t", ctg_num])
                run_ref.log(logstring)
                # set inputs and outputs
                q_contig = ctgs_fas_dir+item
                file_list = (ref_ctg_file, q_contig)
                mauve_outfile = mauve_dir+ctg_num+".mauve"
                aln_segs_dir = aln_segs_root+ctg_num+"/"
                ensure_dir([aln_segs_dir])
                segfile = aln_segs_dir+ctg_num+"_"+ref_n+"_segs.txt"
                open(segfile, 'w').write('')
                # do Mauve alignment
                try:
                    open(ref_ctg_file, 'r')
                    open(q_contig, 'r')
                except IOError:
                    msg = "\nERROR: File missing, cannot align\n\t\t\t"
                    run_ref.log(msg)
                    print msg
                else:
                    align_mauve(file_list, mauve_outfile, mauve_exec)
                    try:
                        # parse Mauve output (without initial clumping)
                        coords = mauver_load2_k0(mauve_outfile+".backbone",
                                                 0, mtype)
                        # chop segments that are too long
                        chop_array = chop_rows(coords, max_size, chop_mode,
                                               mtype)
                        # make detailed pairwise alignments of the segments
                        ref_rec = load_genbank(ref_ctg_file)
                        query_rec = load_fasta(q_contig)
                        iter_align(chop_array, ref_rec, query_rec,
                                   aln_segs_dir, segfile)
                    except IOError:
                        msg = "\nERROR: Mauve alignment failed\n\t\t\t"
                        run_ref.log(msg)
                        print msg
                    except Exception:
                        msg = "\nERROR: Iteration failed\n\t\t\t"
                        run_ref.log(msg)
                        print msg
        print ""

def align_cstrct2ref(run_ref, run_id, timestamp, r_root_dir, run_dirs,
                     genomes, max_size, chop_mode, mtype, mauve_exec):
    """Align constructs pairwise to the reference contig."""
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    ref_ctg_file = run_ref.file
    mauve_root = run_root+run_dirs['mauve_out_dir']+ref_n+"/constructs/"
    segments_root = run_root+run_dirs['aln_seg_dir']+ref_n+"/constructs/"
    scaff_root = run_root+run_dirs['scaffolds_dir']+ref_n+"/"
    ensure_dir([segments_root])
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Align scaffold constructs to reference @",
                         timestamp, "\n"])
    run_ref.log(logstring)
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        scaff_gbk = scaff_root+g_name+"_"+ref_n+"_scaffold.gbk"
        file_list = (ref_ctg_file, scaff_gbk)
        print "\t", g_name, "...",
        # log
        logstring = "".join(["\n", g_name])
        run_ref.log(logstring)
        # set outputs
        mauve_dir = mauve_root+g_name+"/"
        aln_segs_dir = segments_root+g_name+"/"
        ensure_dir([mauve_dir, aln_segs_dir])
        mauve_outfile = mauve_dir+g_name+"_"+ref_n+".mauve"
        segfile = aln_segs_dir+g_name+"_"+ref_n+"_segs.txt"
        # abort if the reference file is not found
        try: open(ref_ctg_file, 'r')
        except IOError:
            msg = "ERROR: Reference file not found"
            print msg
            run_ref.log(msg)
            raise
        # abort if there is no scaffold construct
        try: open(scaff_gbk, 'r')
        except IOError:
            msg = "WARNING: No scaffold construct to align"
            print msg
            run_ref.log(msg)
        else:
            # prep segments file
            open(segfile, 'w').write('')
            # purge any pre-existing sslist file
            sslist_file = scaff_gbk+".sslist"
            if os.path.isfile(sslist_file):
                try: os.remove(sslist_file)
                except Exception: raise
            # do Mauve alignment
            align_mauve(file_list, mauve_outfile, mauve_exec)
            try:
                # parse Mauve output (without initial clumping)
                coords = mauver_load2_k0(mauve_outfile+".backbone", 0, mtype)
                print len(coords), '->',
                logstring = "".join(["\t", str(len(coords))])
                run_ref.log(logstring)
                # chop segments that are too long
                chop_array = chop_rows(coords, max_size, chop_mode, mtype)
                print len(chop_array), 'segments <', max_size, 'bp',
                logstring = "".join(["\t", str(len(chop_array))])
                run_ref.log(logstring)
                # make detailed pairwise alignments of the segments
                ref_rec = load_genbank(ref_ctg_file)
                query_rec = load_genbank(scaff_gbk)
                id = iter_align(chop_array, ref_rec, query_rec,
                                aln_segs_dir, segfile)
                print "@", id, "% id. overall"
                logstring = "".join(["\t", str(id)])
                run_ref.log(logstring)
            except IOError:
                msg = "\nERROR: Mauve alignment failed"
                run_ref.log(msg)
                print msg

def iter_align(coord_array, ref_rec, query_rec, aln_dir, segs_file):
    """Iterate through array of coordinates to make pairwise alignments."""
    # set up the root subdirectories
    seqs = aln_dir+"input_seqs/"
    alns = aln_dir+"output_alns/"
    ensure_dir([seqs, alns])
    aln_id = 0
    aln_len = 0
    # cycle through segments
    for segment_pair in coord_array:
        xa, xb, xc, xd = segment_pair
        # extract the corresponding sequence slices
        ref_seq = ref_rec[abs(xa):abs(xb)]
        query_seq = query_rec[abs(xc):abs(xd)]
        # reverse-complement sequences with negative sign
        if xa < 0 :
            ref_seq = ref_seq.reverse_complement()
        if xc < 0 :
            query_seq = query_seq.reverse_complement()
        # write sequences to file
        mscl_in = seqs+str(xa)+"_"+str(xb)+"_"+str(xc)+"_"+str(xd)+".fas"
        write_fasta(mscl_in, [ref_seq, query_seq])
        # skip segments that are too small to align
        if abs(abs(xa)-abs(xb)) < 10:
            idp = 0
        else:
            # set up outfiles
            mscl_out = alns+str(xa)+"_"+str(xb)+"_"+str(xc)+"_"+str(xd)+".aln"
            logfile = aln_dir+"muscle_log.txt"
            # perform alignment
            align_muscle(mscl_in, mscl_out, logfile)
            idntot = parse_clustal_idstars(mscl_out)
            idp = int((float(idntot)/len(query_seq))*100)
            aln_id += idntot
            aln_len += len(query_seq)
        # write details out to segments file
        line = "\t".join([str(xa), str(xb), str(xc), str(xd), str(idp)+"\n"])
        open(segs_file, 'a').write(line)
    overall_id = int((float(aln_id)/aln_len)*100)
    return overall_id
