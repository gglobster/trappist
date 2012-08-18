import re
import numpy as np
from os import path, listdir
from shutil import copyfile
from loaders import load_genbank, load_multifasta
from writers import write_genbank, write_fasta
from string_ops import multisplit_finder
from common import ensure_dir
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from annotation import annot_ref
from parsing import mauver_load2_k0
from array_tetris import get_anchor_loc
from reporting import ctg_stats
from classes import Reference
from blasting import make_ref_DB

def unpack_genomes(genome, separator, fixed_dirs, ctg_thresholds):
    """Unpack genome files.

    Here, unpacking means extracting data and producing specific files to
    standardize how the information is made available to downstream analysis.
    Depending on the input file format, different unpacking methods are
    invoked. In all cases, this ensures that for each genome, there is a
    multifasta file of the contigs all together as well as a separate Genbank
    file for each contig.

    Supported input file formats are the following:
    - mfas: Basic whole genome sequence in multifasta file of contigs. This
    can be used to process a finished genome in a single Fasta file as well.
    - cgbk: All contigs concatenated in a single GenBank file (Genoscope,
    French WGS). This can be used to process a finished genome in a single
    GanBank file as well.
    # TODO: provide support for other possible input formats

    Unpacking 'cgbk' genomes involves an initial step to detect occurrences
    of the sequence separator and collect the start and stop coordinates of
    each contig. Each pair of coordinates can then be used to extract the
    contig sequence and create a SeqRecord for that contig, which SeqIO
    normally does when it unpacks multifasta files.

    """
    # set up inputs
    infile = genome['file'] #TODO: make GUI input loader (upstream)
    inpath = fixed_dirs['ori_g_dir']+infile
    g_name = genome['name']
    print " ", g_name, "...",
    # prep output destinations
    mfas_dir = fixed_dirs['mfas_contigs_dir']
    fas_dir = fixed_dirs['fas_contigs_dir']+g_name+"/"
    ensure_dir([mfas_dir, fas_dir])
    mfas_file = mfas_dir+g_name+"_contigs.fas"
    records = []
    # select unpacking method
    if genome['input'] is 'fas':
        try: path.exists(inpath) is True
        except ValueError: raise Exception("Bad input file path")
        genome_recs = load_multifasta(inpath)
        # generate GenBank files
        counter = 0
        for rec in genome_recs:
            counter +=1
            ctg_num = str(counter)
            new_id = g_name+"_"+ctg_num  # workaround for long ids
            new_seq = rec.seq
            new_seq.alphabet = generic_dna
            new_rec = SeqRecord(seq=new_seq, id=new_id)
            records.append(new_rec)  # for multifasta output
            fas_file = fas_dir+new_id+".fas"
            write_fasta(fas_file, new_rec)
    elif genome['input'] is 'gbk':
        # load in genome data
        genome_rec = load_genbank(inpath)
        g_string = genome_rec.seq
        # find split coordinates
        coord_pairs = multisplit_finder(g_string, separator)
        # split record
        counter = 0
        for (start, stop) in coord_pairs:
            counter +=1
            ctg_num = str(counter)
            new_record = genome_rec[start:stop]
            new_record.id = g_name+"_"+ctg_num
            records.append(new_record)  # for multifasta output
            fas_file = fas_dir+g_name+"_"+ctg_num+".fas"
            write_fasta(fas_file, new_record)
    else:
        xmsg = "Input file format "+genome['input']+" unspecified/unsupported"
        raise Exception(xmsg)
    print counter, "contigs"
    # write master file
    write_fasta(mfas_file, records)
    # pass records to stats logger
    ctg_stats(g_name, fixed_dirs, ctg_thresholds, records)

def process_ref(ref, ref_annot_flag, r_root_dir, fixed_dirs, run_dirs,
                run_id, timestamp, prot_db_name, project_id):
    """Re-annotate contig and extract reference segments using coordinates."""
    # set inputs and outputs
    run_root = r_root_dir+run_id+"/"
    ref_name = ref['name']
    in_file = fixed_dirs['ori_g_dir']+ref['file']
    seg_out_root = run_root+run_dirs['ref_seg_dir']+ref_name+"/"
    gen_fas_root = fixed_dirs['fas_contigs_dir']+ref_name+"/"
    if ref_annot_flag:
        ref_gbk = run_root+run_dirs['ref_gbk_dir']+ref_name+"_re-annot.gbk"
    else: ## bypass re-annotated ONLY IF ORIGINAL INPUT IS GBK #todo: fix
        ref_gbk = in_file
    ref_fas = run_root+run_dirs['ref_fas_dir']+ref_name+".fas"
    genome_fas = gen_fas_root+ref_name+"_1.fas"
    report_root = run_root+run_dirs['reports']+ref_name+"/"
    ref_log = report_root+run_id+"_"+ref_name+"_log.txt"
    ensure_dir([seg_out_root, report_root, gen_fas_root])
    print " ", ref_name, "...",
    # initialize run_ref object
    run_ref = Reference(ref_name, in_file, ref['input'], ref['seg_mode'],
                        ref['capture'], ref_fas, ref_gbk, seg_out_root,
                        ref_log)
    # initialize reference log
    cl_header = ["# Console log:", run_id, "/", ref_name, timestamp, "\n\n"]
    open(ref_log, 'w').write(" ".join(cl_header))
    # open record and ensure we have a fasta in the right place
    if not path.exists(ref_fas):
        if run_ref.input == 'fas':
            copyfile(in_file, ref_fas)
        elif run_ref.input == 'gbk':
            record = load_genbank(in_file)
            record.id = ref_name
            write_fasta(ref_fas, record)
        else:
            msg = "ERROR: Input not recognized for "+ref_name
            run_ref.log(msg)
            raise Exception(msg)
    # make a BLAST DB
    make_ref_DB(ref, run_id, fixed_dirs, r_root_dir, run_dirs)
    copyfile(ref_fas, genome_fas)
    # re-annotate ref contig
    if ref_annot_flag:
        record = annot_ref(ref_name, ref_fas, prot_db_name, fixed_dirs,
                           project_id)
    else: ## bypass re-annotation ONLY IF ORIGINAL INPUT IS GBK #todo: fix
        record = load_genbank(in_file)
    # load or generate segment definitions
    if run_ref.seg_mode == 'chop':
        run_ref.get_segs_from_chop(len(record.seq), ref['chop_size'])
    elif run_ref.seg_mode == 'list':
        run_ref.get_segs_from_list(ref['segs'])
    elif run_ref.seg_mode == 'feats':
        run_ref.get_segs_from_feats(ref['feat_type'])
    # extract segment sequences
    rec_annot = run_ref.extract_segs_seqs(record, seg_out_root)
    # write re-annotated reference sequence to file
    write_genbank(ref_gbk, rec_annot)
    # report results
    logstring = " ".join([str(len(run_ref.segs)), "segments"])
    print logstring
    run_ref.log(logstring)
    return run_ref

def build_scaffolds(run_ref, r_root_dir, run_dirs, prox_D, separator,
                    genomes, run_id, timestamp, mtype, mode):
    """Build a scaffold of contigs based on the reference.

    This takes contigs that gave positive hits when blasted with reference
    segments. The contigs were aligned against the complete reference in a
    previous step for mapping purposes. Now the output of that step is re-used
    determine their position. A caveat is that if there are natural local
    rearrangements in the sequence relative to the reference, they may not be
    resolved appropriately. The problem is somewhat moderated by the fact that
    this function takes the best (usually the largest) hit region as "anchor"
    to position the contig within the scaffold. But if the rearranged region
    takes up a significant portion of the contig length, the anchoring will
    probably not be called correctly. Visual inspection of the finalized
    maps should help diagnose any such problems. The order can be fixed
    manually using the Mauve Contig Mover, which is part of Mauve 2.

    Note that not all hit contigs are "real" hits, so filtering should be
    applied before scaffolding to generate constructs.

    Model-based filtering produces a list of contigs that will be passed to
    the scaffolder. If filtering manually by looking at the maps,
    there are two options available: either select exclusively OR exclude a
    subset of contigs for the scaffolding process. This is done by listing
    their ID number in the genome dictionaries in the config file then
    resuming the pipeline from this step.

    """
    # set inputs and outputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    ctgs_root = run_root+run_dirs['run_gbk_ctgs_dir']+ref_n+"/"
    mauve_root = run_root+run_dirs['mauve_out_dir']+ref_n+"/contigs/"
    scaffolds_dir = run_root+run_dirs['scaffolds_dir']+ref_n+"/"
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Build scaffold constructs @", timestamp, "\n"])
    run_ref.log(logstring)
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        ctgs_dir = ctgs_root+g_name+"/"
        print "\t", g_name, "...",
        # log
        logstring = "".join(["\n", g_name])
        run_ref.log(logstring)
        # set outputs
        mauve_dir = mauve_root+g_name+"/"
        ensure_dir([mauve_dir, scaffolds_dir])
        scaff_fas = scaffolds_dir+g_name+"_"+ref_n+"_scaffold.fas"
        scaff_gbk = scaffolds_dir+g_name+"_"+ref_n+"_scaffold.gbk"
        # list genbank files in matches directory
        dir_contents = listdir(ctgs_dir)
        anchors_array = np.zeros(1, dtype=[('ctg', 'i4'),
                                           ('start', 'i4'),
                                           ('end', 'i4'),
                                           ('orient', 'i2')])
        # identify contigs we want to select
        subset = []
        for item in dir_contents:
            pattern = re.compile(r'.*_(\d*)\.gbk$')
            match = pattern.match(item)
            if match:
                ctg_num = match.group(1)

                if mode == "exclude":
                    try:
                        if int(ctg_num) in genome[mode]:
                            msg = "("+ctg_num+")"
                            print msg,
                            run_ref.log(msg)
                        else:
                            subset.append(ctg_num)
                    except KeyError:
                        msg = "WARNING: no ignored segments list, including all"
                        print msg
                        msg = ctg_num
                        print msg,
                        subset.append(ctg_num)
                        run_ref.log(msg)
                elif mode == "select":
                    try:
                        if int(ctg_num) in genome[mode]:
                            msg = ctg_num
                            print msg,
                            run_ref.log(msg)
                            subset.append(ctg_num)
                        else:
                            msg = "("+ctg_num+")"
                            print msg,
                            run_ref.log(msg)
                    except KeyError:
                        msg = "WARNING: no selected segments list, including all"
                        print msg
                        msg = ctg_num
                        print msg,
                        subset.append(ctg_num)
                        run_ref.log(msg)
        # at this point we should have a subset of contigs selected
        for ctg_num in subset:
            logstring = "".join(["\t", ctg_num])
            run_ref.log(logstring)
            # set inputs
            mauve_file = mauve_dir+ctg_num+".mauve"
            bb_file = mauve_file+".backbone"
            try:
                # parse Mauve output
                coords = mauver_load2_k0(bb_file, prox_D, mtype)
                # determine which segment to use as anchor
                anchor_seg = get_anchor_loc(coords)
                anchors_array = np.insert(anchors_array, 0,
                                          (ctg_num,
                                           anchor_seg['start'],
                                           anchor_seg['end'],
                                           anchor_seg['orient']))
            except IOError:
                msg = "\tERROR: Mauve alignment not found\n\t"
                print msg
                run_ref.log(msg)
            except Exception:
                msg = "\tERROR: Iteration failure\n\t"
                print msg
                run_ref.log(msg)

        # abort if there is no valid contig to proceed with
        try:
            assert len(anchors_array) > 1 # always 1 left from stub
        except AssertionError:
            msg = "\tWARNING: Contig list empty\n\t"
            print msg
            run_ref.log(msg)
        else:
            # order contigs by anchor location
            anchors_array = np.sort(anchors_array, order='start')
            # load contig records from the genbank files in the matches directory
            ctg_list = []
            for ctg_anchor in anchors_array:
                ctg_num = ctg_anchor['ctg']
                if ctg_num > 0:
                    contig_gbk = ctgs_dir+g_name+"_"+str(ctg_num)+".gbk"
                    record = load_genbank(contig_gbk)
                    if ctg_anchor['orient'] == -1: # flip record
                        record = record.reverse_complement(id=True, name=True,
                            annotations=True, description=True)
                    ctg_list.append(record)
                else: # workaround for having 0 value leftover from stub
                    pass # having it might come in handy in later dev
            # output scaffold files
            write_fasta(scaff_fas, ctg_list)
            scaff_record = SeqRecord('', id='temp')
            scaff_bumper = SeqRecord(separator, id='join')
            for record in ctg_list:
                feat_start = len(scaff_record.seq)
                scaff_record += record
                feat_stop = len(scaff_record.seq)
                scaff_record += scaff_bumper
                feat_loc = FeatureLocation(feat_start, feat_stop)
                pattern = re.compile(r'.*_(\d*)$')
                match = pattern.match(record.id)
                try: ctg_num = match.group(1)
                except Exception: ctg_num = 'N'
                feature = SeqFeature(location=feat_loc,
                                     type='contig',
                                     qualifiers={'id': ctg_num})
                scaff_record.features.append(feature)
            scaff_record.description = g_name+" scaffold from "+ref_n
            try:
                scaff_record.id = g_name
                write_genbank(scaff_gbk, scaff_record[:-100]) # rm last bumper
            except ValueError:
                scaff_record.id = g_name[:10]
                write_genbank(scaff_gbk, scaff_record[:-100]) # rm last bumper
            print ""

def add_refs_2g(genomes, references):
    """ Add references to persistent genomes list."""
    for ref in references:
        genomes.append(ref)
    return genomes