## script to predict orfs on sequences

import re
from os import path
from sys import argv
from libs.common import from_dir, ensure_dir, fas2gbk, gbk2fas, write_genbank, \
    load_genbank, train_prodigal, run_prodigal, load_multifasta
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

origin_dir = "data/"+argv[1]+"/"
file_ext = argv[2]

annot_gbk_dir = origin_dir+"annot_gbk/"
annot_aa_dir = origin_dir+"annot_aa/"
trn_file = origin_dir+"prodigal.trn"

ensure_dir([annot_gbk_dir, annot_aa_dir])

filenames = from_dir(origin_dir, re.compile(r'.*\.'+file_ext+'.*'))

for filename in filenames:
    rec_name = filename[:filename.find("."+file_ext)]

    print rec_name, "...",

    # load data
    if file_ext == 'fas':
        fas_file = origin_dir+"/"+filename
        gbk_file = fas2gbk(fas_file)
        record = load_genbank(gbk_file)
    else:
        gbk_file = origin_dir+"/"+filename
        fas_file = gbk2fas(gbk_file)
        record = load_genbank(gbk_file)

    # run prediction
    annot_aa = annot_aa_dir+rec_name+"_ann.fas"
    annot_gbk = annot_gbk_dir+rec_name+"_ann.gbk"
    if not path.exists(trn_file):
        train_prodigal(fas_file, trn_file, "-q")
    if not path.exists(annot_aa):
        run_prodigal(fas_file, annot_gbk, annot_aa, trn_file, "-q")

    # collect orfs
    record.features = []
    aa_record = load_multifasta(annot_aa)
    counter = 1
    for aa_rec in aa_record:
        this_prot = rec_name+"_"+str(counter)
        # get feature details from description line
        # because prodigal output fails to load as valid genbank
        defline = aa_rec.description
        pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
        match = pattern.match(defline)
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        strand_pos = int(match.group(3))
        feat_loc = FeatureLocation(start_pos, end_pos)
        l_tag = rec_name+"_"+str(counter)
        # consolidation feature annotations
        quals = {'note': defline, 'locus_tag': l_tag,
                 'translation': aa_rec.seq}
        feature = SeqFeature(location=feat_loc,
                             strand=strand_pos,
                             id='cds_'+str(counter),
                             type='CDS',
                             qualifiers=quals)
        record.features.append(feature)
        counter +=1

    # add annotations for Nx100 spacers
    sequence = str(record.seq)
    separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    regex = re.compile(separator, re.IGNORECASE)
    spacers = [match.start() for match in regex.finditer(str(record.seq))]
    #print spacers
    for spacer in spacers:
        space_loc = FeatureLocation(spacer, spacer+100)
        feature = SeqFeature(location=space_loc,
                             type='spacer')
        record.features.append(feature)

    # save record with annotations
    record.description = rec_name+"_with_ORFs"
    record.name = rec_name
    record.dbxrefs = ["Project: "+argv[1]+"/"+rec_name]
    record.seq.alphabet = generic_dna
    write_genbank(annot_gbk, record)

    print "OK"


