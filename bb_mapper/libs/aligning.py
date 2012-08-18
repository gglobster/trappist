import os, re
import subprocess
from os import path
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from writers import ensure_dir, write_fasta
from parsing import parse_clustal_idstars, mauver_load2_k0
from loaders import load_genbank, from_dir
from array_tetris import chop_rows

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

def mauve_pw_align(ref, query, dirs, run, max_size, chop_mode, mauve_exec):
    """Set up and perform a pairwise alignment with Mauve."""
    aln_dir = dirs['root']+run+dirs['aln_segs']
    mauve_dir = dirs['root']+run+dirs['mauve']
    # set outputs
    mauve_outfile = mauve_dir+ref.name+"_"+query.name+".mauve"
    segfile = aln_dir+ref.name+"_"+query.name+"_segs.txt"
    # check for existing alignment
    if path.exists(segfile):
        print "already done"
    else:
        # prep segments file
        open(segfile, 'w').write('')
        # purge any pre-existing sslist files
        sslist_files = from_dir(dirs['seqfiles'], re.compile(r'.*\.sslist.*'))
        for sslist in sslist_files:
            try: os.remove(dirs['seqfiles']+sslist)
            except Exception: raise
        # do Mauve alignment
        file_list = [ref.gbk, query.gbk]
        align_mauve(file_list, mauve_outfile, mauve_exec)
        try:
            # parse Mauve output (without initial clumping)
            coords = mauver_load2_k0(mauve_outfile+".backbone", 0)
            print "\nSegment results:", len(coords), '->',
            # chop segments that are too long
            chop_array = chop_rows(coords, max_size, chop_mode)
            print len(chop_array), 'segments <', max_size, 'bp'
            # make detailed pairwise alignments of the segments
            print "Aligning segments ..."
            ref_rec = load_genbank(ref.gbk)
            query_rec = load_genbank(query.gbk)
            id = iter_align(chop_array, ref_rec, query_rec,
                            dirs['aln_segs'], segfile)
            print "Results:", id, "% id. overall"
        except IOError:
            print "\nERROR: Mauve alignment failed"

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
        print ".",
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
    print ""
    return overall_id

