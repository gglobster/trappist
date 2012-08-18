## script to predict orfs on sequences

import re
from os import path
from sys import argv
from libs.common import from_dir, ensure_dir, fas2gbk, gbk2fas, write_genbank, \
    load_genbank, train_prodigal, run_prodigal, load_multifasta,\
    remote_blastp_2file, collect_topNhits
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

origin_dir = "data/"+argv[1]+"/"
seq_dir = origin_dir+argv[2]+"/"
file_ext = argv[3]

if len(argv) < 5:
    trim_ids = ''
else:
    trim_ids = argv[4]

blast_dir = origin_dir+"blast/"
hits_dir = origin_dir+"hits/"
remote_prot_db = "nr"

annot_gbk_dir = origin_dir+"annot_gbk/"
annot_aa_dir = origin_dir+"annot_aa/"
trn_file = origin_dir+"prodigal.trn"

ensure_dir([annot_gbk_dir, annot_aa_dir, blast_dir, hits_dir])

filenames = from_dir(seq_dir, re.compile(r'.*\.'+file_ext+'.*'))

for filename in filenames:
    rec_name = filename[:filename.find(trim_ids+"."+file_ext)]

    print rec_name, "..."

    # load data
    if file_ext == 'fas':
        fas_file = seq_dir+"/"+filename
        gbk_file = fas2gbk(fas_file)
        record = load_genbank(gbk_file)
    else:
        gbk_file = seq_dir+"/"+filename
        fas_file = gbk2fas(gbk_file)
        record = load_genbank(gbk_file)

    assert record.id

    # run prediction
    annot_aa = annot_aa_dir+rec_name+"_ann.fas"
    annot_gbk = annot_gbk_dir+rec_name+"_ann.gbk"
    if not path.exists(trn_file):
        train_prodigal(fas_file, trn_file, "-q")
    if not path.exists(annot_aa):
        run_prodigal(fas_file, annot_gbk, annot_aa, trn_file, "-q")

    # blast the protein sequences against the remote DB
    record.features = []
    evalue = 0.01
    proteins = load_multifasta(annot_aa)
    for protein in proteins:
        print "  ", protein.id
        rec_hits_dir = hits_dir+rec_name+"/"
        ensure_dir([rec_hits_dir])
        hits_out = open(rec_hits_dir+protein.id+".txt", 'w')
        hits_out.write(" ".join([protein.id, "vs.", remote_prot_db,
                                 "@evalue =", str(evalue), "\n"]))
        temp_out = remote_blastp_2file(protein.seq, remote_prot_db,
                                       blast_dir+rec_name+"_temp.xml",
                                       evalue)
        #temp_out = blast_dir+rec_name+"_temp.xml"
        # collect best 10 hits
        rec_hits = collect_topNhits(temp_out, 10)
        for hit in rec_hits:
            if hasattr(hit, 'hsps'):
                hit_cleaned = re.sub(r">", "\n", str(hit))
                hit_cleaned = re.sub(r"Length = \d+", "", hit_cleaned)
                hits_out.write("\n\n"+hit_cleaned)
                for hsp in hit.hsps:
                    hits_out.write(str(hsp)+"\n")
            else:
                hits_out.write("\nNo hits\n")
                break
        hits_out.close()
        # get predicted feature details from description line
        # because prodigal output fails to load as valid genbank
        defline = protein.description
        pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
        match = pattern.match(defline)
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        strand_pos = int(match.group(3))
        feat_loc = FeatureLocation(start_pos-1, end_pos) # adjust for 0-index
        l_tag = protein.id
        # consolidate feature annotations
        quals = {'note': defline, 'locus_tag': l_tag,
                 'translation': protein.seq}
        feature = SeqFeature(location=feat_loc,
                             strand=strand_pos,
                             id=protein.id,
                             type='CDS',
                             qualifiers=quals)
        record.features.append(feature)

    # save record with annotations
    record.description = rec_name+"_with_ORFs"
    record.name = rec_name
    record.dbxrefs = ["Project: "+argv[1]+"/"+rec_name]
    record.seq.alphabet = generic_dna
    write_genbank(annot_gbk, record)

    print "OK"


