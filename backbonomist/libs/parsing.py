import re, numpy as np
from os import listdir, path
from shutil import copyfile
from loaders import read_array, td_txt_file_load, load_fasta
from writers import write_fasta
from common import ensure_dir
from array_tetris import extract_nonzero, clump_rows
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

def glompX_blast_out(genomes, run_ref, blast_mode, r_root_dir, run_dirs,
                     run_id, fixed_dirs, blast_dtypes, references,
                     min_nt_match, min_nt_score, min_nt_idp, min_aa_match,
                     min_aa_score, min_aa_idp, capture_span, timestamp):
    """Collect Blast results and extract match contigs."""
    # load inputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    match_root = run_root+run_dirs['match_out_dir']+ref_n+"/"
    capture_root = run_root+run_dirs['capture_dir']+ref_n+"/"
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Collect Blast results @", timestamp, "\n\n"])
    run_ref.log(logstring)
    # collect results
    ref_hits = {}
    control_scores = []
    run_ref.log("Segs/Gs\t")
    run_ref.log("\t".join([genome['name'] for genome in genomes]))
    for seg in run_ref.segs:
        seg_n = seg['name']
        print "\t", seg_n, "...",
        run_ref.log("".join(["\n", seg_n]))
        blast_dir = run_root+run_dirs['blast_out_dir']+ref_n+"/"+seg_n+"/"
        capture_dir = capture_root+"/"+seg_n+"/"
        ensure_dir([blast_dir, capture_dir])
        ref_flag = True
        for genome in genomes:
            g_name = genome['name']
            print "|",
            # process
            if g_name not in ref_hits.keys():
                ref_hits[g_name] = {}
            matches_dir = match_root+g_name+"/"
            ensure_dir([matches_dir])
            blast_infile = blast_dir+g_name+"_out.txt"
            genome_ctg_dir = fixed_dirs['fas_contigs_dir']+g_name+"/"
            try:
            	rec_array = read_array(blast_infile, blast_dtypes)
            except Exception:
            	rec_array = []
            if len(rec_array) > 0:  # take qualified hits
                p_cnt = 0
                n_cnt = 0
                if g_name in [ref['name'] for ref in references]:
                    copyfile(genome_ctg_dir+g_name+"_1.fas",
                             matches_dir+g_name+".fas")
                    if ref_flag:
                        # positive control TODO: better solution
                        control_scores.append(rec_array[0][11])
                        ref_flag = False
                for line in rec_array:
                    idp = line[2]
                    q_start, q_stop = line[8], line[9]
                    score = line[11]
                    length = abs(q_stop-q_start)
                    # check the blast mode to use the right thresholds
                    if blast_mode == 'n' or blast_mode == 'tx':
                        min_match = min_nt_match
                        min_score = min_nt_score
                        min_idp = min_nt_idp
                    elif blast_mode == 'tn':
                        min_match = min_aa_match
                        min_score = min_aa_score
                        min_idp = min_aa_idp
                    else: # default to nucleotide mode
                        min_match = min_nt_match
                        min_score = min_nt_score
                        min_idp = min_nt_idp
                    if length>min_match and score>min_score and idp>min_idp:
                        print "+",
                        p_cnt +=1
                        contig_id = line[1]
                        if contig_id not in ref_hits[g_name].keys():
                            ref_hits[g_name][contig_id] = {seg_n: score}
                        else:
                            ref_hits[g_name][contig_id][seg_n] = score
                        pattern = re.compile(r'('+contig_id+')\.fas')
                        for item in listdir(genome_ctg_dir):
                            match = re.match(pattern, item)
                            if match:
                                fas_file = matches_dir+match.group(1)+".fas"
                                if not path.exists(fas_file):
                                    copyfile(genome_ctg_dir+item, fas_file)
                        # context capture
                        capture_flag = False
                        while True:
                            try:
                                if int(seg_n) in run_ref.capture:
                                    capture_flag = True
                                else:
                                    break
                            except ValueError:
                                if seg_n in run_ref.capture:
                                    capture_flag = True
                                else:
                                    break
                            else:
                                break
                        if capture_flag:
                            # load the sequence
                            contig_file = matches_dir+contig_id+".fas"
                            contig_rec = load_fasta(contig_file)
                            # check orientation
                            if q_start < q_stop:
                                c_start = q_start-capture_span
                                c_stop = q_stop+capture_span
                            else:
                                c_start = q_stop-capture_span
                                c_stop = q_start+capture_span
                            print c_start, c_stop
                            # check limits
                            if c_start < 0:
                                c_start = 1
                            if c_stop > len(contig_rec.seq):
                                c_stop = len(contig_rec.seq)
                            # proceed
                            cxt_file = capture_dir+g_name+"_"+contig_id+".fas"
                            cxt_rec = SeqRecord(id=contig_id+"_"
                                                    +str(c_start)+"_"
                                                    +str(c_stop),
                                                seq=contig_rec.seq
                                                    [c_start:c_stop])
                            write_fasta(cxt_file, cxt_rec)
                    else:
                        print "-",
                        n_cnt +=1
                if n_cnt > 0:
                    logstring = "".join(["\t", str(p_cnt), " (",
                                         str(n_cnt), ")"])
                else:
                    logstring = "".join(["\t", str(p_cnt)])
                run_ref.log(logstring)
            else:
                print "-",
                run_ref.log("".join(["\t", "0"]))
        print ""
    return ref_hits, control_scores

def mauver_load2_k0(file, threshold, mtype):
    """Parse Mauve coordinates file to extract segment coordinates.

    This loads the coordinates data into a Numpy array. All rows that contain
    a zero value are deleted from the array. Rows are then collapsed if the
    coordinates are close (under a user-specified threshold) and in the same
    orientation. The resulting array is then returned.

    Important note: this function only works for pairwise alignments!

    """
    # load file data into numpy array
    raw_array = np.loadtxt(file, skiprows=1, dtype=mtype)
    # set up an array stub to receive non-zero rows (leaves a (1,1) pair)
    stub_array = np.ones(1, dtype=mtype)
    # eliminate rows containing elements of value 0
    try:
        nz_array = extract_nonzero(raw_array, stub_array)
    except TypeError:
        nz_array = np.append(stub_array, raw_array)
    # collapse rows
    cl_array = clump_rows(nz_array, threshold, mtype)
    return cl_array

def collect_cogs(blast_out):
    """Collect hits from Blast XML (not just for COGs anymore)."""
    results = {}
    blast_records = NCBIXML.parse(open(blast_out))
    for record in blast_records:
        rec_key = record.query_id
        results[rec_key+'_def'] = record.query
        if record.alignments: # ignores searches with no hits
            top_hit = record.alignments[0]
            results[rec_key] = top_hit.title
        else:
            results[rec_key] = 'no match'
    return results

def parse_clustal_idstars(filename):
    """Parse ClustalW output file to estimate identity percentage."""
    #from analysis.text_manipulation import td_txt_file_load
    raw_lines = td_txt_file_load(filename, 3)
    single_line = ''.join(raw_lines)
    idnstar = re.compile(r"\*")
    idnalls = re.findall(idnstar, single_line)
    idntot = len(idnalls)
    return idntot # total number of identical nucleotide positions
