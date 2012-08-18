from os import path
from common import ensure_dir
from plotting import plot_ctg_stats, hits_heatmap_multi
import numpy as np

def log_start_run(run_id, base_root, project_id, run_dirs, timestamp):
    """Record run initiated in the dataset log."""
    set_log = base_root+"/"+project_id+"_log.html"
    param_file = run_id+"/"+run_dirs['reports']+run_id+"_dataset.txt"
    header = "<p><b>Log of processing runs for "+project_id+"</b></p><p><ul>"
    set_log_htm = ["<li><b>", run_id, "</b>&nbsp;- initiated ", timestamp,
                   " (<a href='", param_file, "'>dataset</a>)</li>"]
    linkline = "".join(set_log_htm)
    if not path.isfile(set_log): # first time running this dataset?
        open(set_log, 'w').write(header)
    open(set_log, 'a').write(linkline)

def log_end_run(run_id, base_root, project_id, timestamp):
    """Record run in the dataset log."""
    set_log = base_root+"/"+project_id+"_log.html"
    run_report = run_id+"/"+run_id+"_report.html"
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- completed ", timestamp,
                   " (<a href='", run_report, "'>report</a>)</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def log_resume_run(run_id, base_root, project_id, timestamp, step):
    """Record run resumed in the dataset log."""
    set_log = base_root+"/"+project_id+"_log.html"
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- resumed ", timestamp,
                   " (from step "+str(step)+")</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def save_datasumm(run_id, blast_mode, r_root_dir, run_dirs, genomes,
                  references, project_id, project_date, timestamp):
    """Save a summary of the dataset composition to file."""
    print " ", run_id,
    report_root = r_root_dir+run_id+"/"+run_dirs['reports']
    param_file = report_root+run_id+"_dataset.txt"
    # genome data
    g_title = "\t".join(["## Genomes"])
    g_header = "\t".join(["## Name", "Offset", "Format", "File"])
    g_data = [g_title, g_header]
    for genome in genomes:
        g_data.append("\t".join([genome['name'], str(genome['offset'][1]),
                                 genome['input'], genome['file']]))
    g_block = "\n".join(g_data)
    # references data
    r_title = "\t".join(["## References"])
    r_header = "\t".join(["## Name", "Segments", "File"])
    r_data = [r_title, r_header]
    for ref in references:
        if ref['seg_mode'] == 'list':
            r_data.append("\t".join([ref['name'], str(len(ref['segs'])),
                                     ref['file']]))
            r_data.append("\t".join(["###", "Segments"]))
            r_data.append("\t".join(["###", "Name", "Start", "Stop"]))
            for seg in ref['segs']:
                r_data.append("\t".join(["", seg['name'],
                                         str(seg['coords'][0]),
                                         str(seg['coords'][1])]))
        else:
            r_data.append("\t".join([ref['name'], ref['file']]))
            if ref['seg_mode'] == 'chop':
                r_data.append("\t".join(["###", "Segments generated from size",
                                     str(ref['chop_size'])]))
            elif  ref['seg_mode'] == 'feats':
                r_data.append("\t".join(["###", "Segments generated from ",
                                         ref['feat_type'], "features"]))
            else:
                r_data.append("\t".join(["###", "Unknown segmenting scheme"]))
    r_block = "\n".join(r_data)
    # text block
    txt = ["# Project", project_id,
           "# Run ID", run_id,
           "# Date generated", project_date,
           "# Date processing initiated", timestamp,
           "# Dataset composition", g_block, r_block,
           "# Blast mode", blast_mode]
    # write to file
    open(param_file, 'w').write("\n".join(txt))
    print "dataset summary saved to file"

def init_reports(run_id, fixed_dirs, ctg_thresholds, timestamp):
    """Record run info in the logs."""
    # set inputs and outputs
    ctg_stats_file = fixed_dirs['ctg_stats']+"contig_stats.txt"
    # initialize contig stats report
    kb1, kb2, kb3 = ctg_thresholds
    cs_header = "# Contig size distribution statistics\n"
    if not path.isfile(ctg_stats_file): # first time running?
        open(ctg_stats_file, 'w').write(cs_header)
    cs_txt = ["#\nRun ", run_id, " initiated ", timestamp,
              "\nGenome", "\t[<", str(kb1), "]",
              "\t[", str(kb1), "-", str(kb2), "]",
              "\t[", str(kb2), "-", str(kb3), "]",
              "\t[>", str(kb3), "]"]
    open(ctg_stats_file, 'a').write("".join(cs_txt))

def ctg_stats(g_name, fixed_dirs, ctg_thresholds, records):
    """Report on size distribution of contigs in a genome."""
    i = 1000
    kb1, kb2, kb3 = ctg_thresholds
    # separate length categories
    small_ctgs = [len(rec)/i for rec in records if len(rec) <kb1*i]
    mid_ctgs = [len(rec)/i for rec in records if kb1*i <= len(rec) < kb2*i]
    big_ctgs = [len(rec)/i for rec in records if kb2*i <= len(rec) <= kb3*i]
    oversized = [len(rec)/i for rec in records if len(rec) >=kb3*i]
    # write to report file
    ctg_stats_file = fixed_dirs['ctg_stats']+"contig_stats.txt"
    txt = ["\n"+g_name, str(len(small_ctgs)), str(len(mid_ctgs)),
           str(len(big_ctgs)), str(len(oversized))]
    open(ctg_stats_file, 'a').write("\t".join(txt))
    # make figure
    fname = fixed_dirs['ctg_stats']+g_name+"_ctg_stats.png"
    ctg_cats = small_ctgs, mid_ctgs, big_ctgs
    plot_ctg_stats(ctg_cats, fname, ctg_thresholds)

def matches_table(match_dict, r_root_dir, run_dirs, timestamp):
    """Compile tables and graphs of contig matches."""
    # unpack results
    run_ref = match_dict['ref']
    run_id = match_dict['run']
    ref_hits = match_dict['hits']
    ctl_scores = match_dict['ctl']
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    report_root = run_root+run_dirs['reports']+ref_n+"/"
    hits_table_txt = report_root+run_id+"_"+ref_n+"_hits_table.txt"
    red_hits_fig = report_root+run_id+"_"+ref_n+"_sum_hits.pdf"
    mf_hits_fig = report_root+run_id+"_"+ref_n+"_all_hits.pdf"
    ensure_dir([report_root])
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Make matches results table & graphs @",
                         timestamp, "\n"])
    run_ref.log(logstring)
    # process
    rep_fhandle = open(hits_table_txt, 'w')
    rep_fhandle.write("# Matches to the reference segments of "+ref_n)
    segs = [seg['name'] for seg in run_ref.segs]
    red_g_list = []
    mf_g_list = []
    mf_ctgs = []
    g_names = []
    # traverse dict of results
    for g_name in sorted(ref_hits, reverse=True):
        print "\t", g_name,
        run_ref.log("".join(["\n", g_name, "\t"]))
        genome_root = report_root+g_name+"/"
        ensure_dir([genome_root])
        g_table_fig = genome_root+g_name+"_hits_"+ref_n+".pdf"
        # collect scores
        g_list = []
        g_header = ["\n\n"+g_name]+segs
        g_lines = ["\t".join(g_header)]
        for contig in ref_hits[g_name]:
            ctg_hits = [ref_hits[g_name][contig][seg]
                        if seg in ref_hits[g_name][contig]
                        else 0 for seg in segs]
            g_list.append(ctg_hits)
            ctg_line = "\t".join([contig]+[str(score) for score in ctg_hits])
            g_lines.append(ctg_line)
        rep_fhandle.write("\n".join(g_lines))
        # normalize scores
        g_array = np.array(g_list)
        if not len(g_list):
            g_norm = np.zeros(len(ctl_scores))
            g_norm = np.reshape(g_norm, (1, len(ctl_scores)))
        else:
            try:
                g_norm = np.divide(g_array, ctl_scores)
            except ValueError:
                msg = "ERROR: problem processing hits for "+g_name
                run_ref.log(msg)
                print msg,
                g_norm = np.zeros(len(ctl_scores))
                g_norm = np.reshape(g_norm, (1, len(ctl_scores)))
        # graph detailed scores per genome
        contigs = ref_hits[g_name].keys()
        hits_heatmap_multi(ref_n, segs, [g_name], [contigs], [g_norm],
                           g_table_fig)
        # prep for global graphs
        if len(g_array) > 1:
            try:
                g_reduce = np.amax(g_array, 0)
                g_reduce = np.divide(g_reduce, ctl_scores)
            except Exception:
                msg = "ERROR: #1 problem preparing global graphs for "+g_name
                run_ref.log(msg)
                print msg,
        else:
            try:
                g_reduce = np.reshape(g_norm, (len(ctl_scores),))
            except Exception:
                msg = "ERROR: #2 problem preparing global graphs for "+g_name
                run_ref.log(msg)
                print msg,
        red_g_list.append(g_reduce)
        mf_g_list.append(g_norm)
        mf_ctgs.append(contigs)
        g_names.append(g_name)
        # estimate total score per genome
        g_score = np.sum(g_reduce)
        print g_score
        run_ref.log(str(g_score))
    # graph detailed scores for all genomes
    hits_heatmap_multi(ref_n, segs, g_names, mf_ctgs, mf_g_list, mf_hits_fig)
    # graph summarized scores for all genomes
    red_g_array = np.array(red_g_list)
    g_names.reverse()
    hits_heatmap_multi(ref_n, segs, [1], [g_names], [red_g_array],
                       red_hits_fig)
    # done
    rep_fhandle.close()

