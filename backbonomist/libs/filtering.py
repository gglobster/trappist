from pprint import pprint
from common import ensure_dir
from loaders import load_fasta
from scipy.stats import nanmean
import numpy as np


def filter_contigs(run_ref, run_id, genomes, norm_matches, seg_size, threshold, r_root_dir, run_dirs, fixed_dirs, timestamp):
    """Filter contigs."""
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    fas_root = fixed_dirs['fas_contigs_dir']
    report_root = run_root+run_dirs['reports']+ref_n+"/"
    ensure_dir([report_root])
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Filter contigs @", timestamp, "\n"])
    run_ref.log(logstring)
    # process
    
    # evaluate segment specificity using negative controls
    neg_controls = [genome['name'] for genome in genomes if ('ctrl' in genome.keys() and genome['ctrl'] == 'neg')]
    neg_dat = [norm_matches[g_name]['ctg_scores'] for g_name in neg_controls]
    neg_RA = np.vstack(neg_dat)
    neg_mean = nanmean(neg_RA, axis=0)
    # process the genomes we're testing
    test_genomes = [genome['name'] for genome in genomes if not ('ctrl' in genome.keys())]
    for g_name in test_genomes: 
        print "\t", g_name,
        ctg_hits = norm_matches[g_name]['ctg_scores']
        ctg_stats = {}
        #process individual contigs
        counter = 0
        for ctg_RA in ctg_hits:
            # identify this contig by name
            ctg_name = norm_matches[g_name]['ctg_names'][counter]
            counter += 1
            # subtract background signal from match scores
            recal_ctg_RA = np.subtract(ctg_RA, neg_mean)
            recal_ctg_RA = recal_ctg_RA.clip(min=0)
            # compute total similarity score
            s_score = np.sum(recal_ctg_RA)
            # compute clustering score (primitive)
            streak = False
            c_score = 0
            for hit in recal_ctg_RA:
                if hit == 0:
                    if streak == True:
                        c_score += -1
                        streak = False
                    else: 
                        c_score += 0
                elif hit > 0:
                    if streak == True:
                        c_score += 2
                    else: 
                        c_score += 1
                        streak = True
            # compute backbone vs. cargo burden
            ctg_rec = load_fasta(fas_root+g_name+"/"+ctg_name+".fas")
            bbone = np.sum(np.ma.make_mask(recal_ctg_RA))*seg_size
            if bbone > len(ctg_rec):
                bbone = len(ctg_rec)    # workaround for last segment being always a little short
            cargo = len(ctg_rec) - bbone
            # make inverted array mask (used for redundancy detection)
            ctg_mask = np.ma.getmaskarray(np.ma.masked_equal(recal_ctg_RA,0))
            # consolidate contig information
            ctg_stats[ctg_name] = {'s_score': s_score, 
                                    'c_score': c_score, 
                                    'vector': recal_ctg_RA, 
                                    'inv_mask':ctg_mask,
                                    'bbone': bbone,
                                    'cargo': cargo}
        # detect redundant contigs
        ### use np.ma.mask_or(m1, m2)
        ### if any elements returns false there is a redundancy between two contigs
        ### if so evaluate which has better c_score and s_score
        
        # compute overall stats for the genome
        gs_score = sum([ctg_stats[contig]['s_score'] for contig in ctg_stats])
        gc_score = sum([ctg_stats[contig]['c_score'] for contig in ctg_stats])
        g_bbone = sum([ctg_stats[contig]['bbone'] for contig in ctg_stats])
        g_cargo = sum([ctg_stats[contig]['cargo'] for contig in ctg_stats])
        print gs_score, gc_score, g_bbone, g_cargo,
        # 
        if gs_score > threshold:
            ## run plotters again 
            ## pass the genome on to the next step (others will be dropped)
            print "MATCH"
        else:
            print "(-)"
        
        