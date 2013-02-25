import os, sys, time
from pprint import pprint

import cPickle as pickle
from datetime import datetime
from libs.common import ensure_dir
from libs.genome_tetris import process_ref, unpack_genomes, add_refs_2g, \
    build_scaffolds
from libs.blasting import make_genome_DB, basic_batch_blast
from libs.parsing import glompX_blast_out
from libs.annotation import annot_genome_contigs
from libs.mapping import prep_maps
from libs.aligning import align_cstrct2ref, align_ctg2ref
from libs.reporting import save_datasumm, log_start_run, log_end_run, \
    init_reports, log_resume_run, matches_table
from libs.filtering import filter_contigs

from run_config import *
from run_sets import references, genomes

i_msg = "\n", \
        "##################################################\n", \
        "### Backbonomist v. 0.9                        ###\n", \
        "### Copyright 2012 Geraldine A. Van der Auwera ###\n", \
        "##################################################\n"

def usage():
    print "Usage: python gobble.py bb_fish -args (don't call this directly!)"

def info():
    print i_msg

def get_arg(args, key):
    return args[args.index(key)+1]

def main(args=None):

    info()

    if args is None:
        args = sys.argv[1:]

    if "-h" in args or "--help" in args:
        usage()
        sys.exit(2)

    if "-id" in args:
        run_id = get_arg(args, "-id")
    else:
        run_id = str(int(time.time())) # use timestamp as unique run identifier

    if "-bm" in args:
        blast_mode = get_arg(args, "-bm")
    else:
        blast_mode = 'n' # nucleotide blast by default

    if "-resume" in args:
        step = int(get_arg(args, "-resume"))
        resume = True
    else:
        step = 0
        resume = False

    if "-short" in args:
        limit = 6
    else:
        limit = 100 # unnecessarily high cap
        
    if "-filter" in args:
        threshold = int(get_arg(args, "-filter"))
        resume = True
        step = 6
        limit = 7
    else:
        threshold = 5 # reduce for small references
        limit = 100 # unnecessarily high cap
        
    if "-ctg" in args:
        ctg_subset = get_arg(args, "-ctg")
    else:
        ctg_subset = 'exclude'

    if "-g" in args:
        g_select = get_arg(args, "-g")
    else:
        g_select = None

    start_timestamp = str(datetime.now())

    # ensure existence of all directories
    ensure_dir(fixed_dirs.values())
    run_dirs_go = ["".join([r_root_dir, run_id, "/", rdir])
                   for rdir in run_dirs.values()]
    ensure_dir(run_dirs_go)

    # define pickle paths
    pickle_root = r_root_dir+run_id+"/"+run_dirs['pickles']+run_id
    ref_pickles = pickle_root+"_refs.p"
    genome_pickles = pickle_root+"_genomes.p"
    blast_pickles = pickle_root+"_blast.p"
    match_pickles = pickle_root+"_matches.p"
    norm_pickles = pickle_root+"_norm.p"

    # check for pickles
    run_refs = []
    run_gs = []
    run_blast = False
    run_matches = []
    run_norm_matches = {}

    if resume:
        # normalized matches
        if step > 5:
            try: run_norm_matches = pickle.load(open(norm_pickles, 'rb'))
            except IOError:
                print "WARNING: Could not load norm pickle"
                run_norm_matches = {}
                step = 5
        # matches
        if step > 4:
            try: run_matches = pickle.load(open(match_pickles, 'rb'))
            except IOError:
                print "WARNING: Could not load matches pickle"
                run_matches = []
                step = 4
        # blast
        if step > 3:
            try: run_blast = pickle.load(open(blast_pickles, 'rb'))
            except IOError:
                print "WARNING: Could not load blast pickle"
                run_blast = False
                step = 3
        # genomes
        if step > 2:
            try: run_gs = pickle.load(open(genome_pickles, 'rb'))
            except IOError:
                print "WARNING: Could not load genomes pickle"
                run_gs = []
                step = 2
        # references
        if step > 1:
            try: run_refs = pickle.load(open(ref_pickles, 'rb'))
            except IOError:
                print "WARNING: Could not load refs pickle"
                run_refs = []
                step = 1
    else:
        step = 0     

    ## pipeline

    print "starting pipeline"

    print step, limit

    if resume:
        log_resume_run(run_id, base_root, project_id, start_timestamp, step)

    else:
        print "\n###", step, ". Set up logging & reporting ###\n"
        log_start_run(run_id, base_root, project_id, run_dirs, start_timestamp)
        save_datasumm(run_id, blast_mode, r_root_dir, run_dirs, genomes,
                      references, project_id, project_date, start_timestamp)
        init_reports(run_id, fixed_dirs, ctg_thresholds, start_timestamp)
        step +=1

    while step < limit:

        if step is 1:
            print "\n###", step, ". Prepare references ###\n"
            for ref in references:
                timestamp = str(datetime.now())
                ref_obj = process_ref(ref, ref_annot_flag, r_root_dir,
                                      fixed_dirs, run_dirs, run_id,
                                      timestamp, prot_db_name, project_id)
                run_refs.append(ref_obj)
            if os.path.exists(ref_pickles):
                os.remove(ref_pickles)
            pickle.dump(run_refs, open(ref_pickles, 'wb'))
            step +=1

        elif step is 2:
            print "\n###", step, ". Prepare genomes ###\n"
            for genome in genomes:
                unpack_genomes(genome, separator, fixed_dirs, ctg_thresholds)
                make_genome_DB(genome, fixed_dirs)
            run_gs = add_refs_2g(genomes, references)
            if os.path.exists(genome_pickles):
                os.remove(genome_pickles)
            pickle.dump(run_gs, open(genome_pickles, 'wb'))
            step +=1

        elif step is 3:
            print "\n###", step, ". Blast reference segments against genomes ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                run_blast = basic_batch_blast(run_gs, ref, blast_mode,
                                              r_root_dir, run_dirs,
                                              fixed_dirs, blast_prefs,
                                              run_id, timestamp)  
                if os.path.exists(blast_pickles):
                   os.remove(blast_pickles)
            pickle.dump(run_blast, open(blast_pickles, 'wb'))
            step +=1

        elif step is 4:
            print "\n###", step, ". Collect Blast results ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                ref_hits, ctl_scores = glompX_blast_out(run_gs, ref,
                                                        blast_mode,
                                                        r_root_dir, run_dirs,
                                                        run_id, fixed_dirs,
                                                        blast_dtypes,
                                                        references,
                                                        min_nt_match, 
                                                        min_nt_score,
                                                        min_nt_idp,
                                                        min_aa_match,
                                                        min_aa_score,
                                                        min_aa_idp,
                                                        capture_span,
                                                        timestamp)
                ref_matches = {'ref': ref, 'run': run_id, 'hits': ref_hits,
                               'ctl': ctl_scores}
                run_matches.append(ref_matches)
            if os.path.exists(match_pickles):
                os.remove(match_pickles)
            pickle.dump(run_matches, open(match_pickles, 'wb'))
            step +=1

        elif step is 5:
            print "\n###", step, ". Make match results table & graphs ###\n"
            for ref_matches in run_matches:
                timestamp = str(datetime.now())
                ref_norm_matches = matches_table(ref_matches, r_root_dir, run_dirs, timestamp)
                run_norm_matches[ref_matches['ref'].name] = ref_norm_matches
            if os.path.exists(norm_pickles):
                os.remove(norm_pickles)
            pickle.dump(run_norm_matches, open(norm_pickles, 'wb'))
            step +=1

        elif step is 6:
            print "\n###", step, ". Filter matching contigs ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                filter_contigs(ref, run_id, genomes, run_norm_matches[ref.name], chop_size, threshold, r_root_dir, 
                            run_dirs, fixed_dirs, timestamp)
            step +=1

        elif step is 7:
            print "\n###", step, ". Annotate matching contigs ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                annot_genome_contigs(ref, prot_db_name, fixed_dirs, r_root_dir,
                                 run_id, run_dirs, genomes, project_id,
                                 timestamp, blast_prefs)
            step +=1

        elif step is 8:
            print "\n###", step, ". Align contigs pairwise to reference ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                align_ctg2ref(ref, run_id, timestamp, r_root_dir, run_dirs, 
                              genomes, mauve_exec, max_size, chop_mode, mtype)
            step +=1

        elif step is 9:
            print "\n###", step, ". Construct backbone-based scaffolds ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                build_scaffolds(ref, r_root_dir, run_dirs, prox_D,
                                separator, genomes, run_id, timestamp,
                                mtype, ctg_subset)
            step +=1

        elif step is 10:
            print "\n###", step, ". Align constructs pairwise to reference ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                align_cstrct2ref(ref, run_id, timestamp, r_root_dir, run_dirs,
                     genomes, max_size, chop_mode, mtype, mauve_exec)
            step +=1

        elif step is 11:
            print "\n###", step, ". Generate maps ###\n"
            for ref in run_refs:
                timestamp = str(datetime.now())
                prep_maps(ref, run_id, timestamp, g_select, r_root_dir,
                          run_dirs, genomes, fixed_dirs, segtype, min_size,
                          fct_flags, fct_colors, idpt)
            step +=1

        elif step > 12:
            break

    stop_timestamp = str(datetime.now())
    log_end_run(run_id, base_root, project_id, stop_timestamp)
    print "\n### Nothing more to do! ###\n"


if __name__ == "__main__":
    sys.exit(usage())