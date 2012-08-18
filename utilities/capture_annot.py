# script to capture features based on annotation tags

import re
from sys import argv
from libs.common import load_genbank, write_fasta, ensure_dir, from_dir

data_dir = "data/"+argv[1]+"/"
dir_in = "data/"+argv[2]+"/"
feat_type = argv[3]
feat_tag = argv[4]
feat_name = argv[5]

main_out = data_dir+feat_name+"_seqs.fas"

records = []
ensure_dir([data_dir])

filenames = from_dir(dir_in, re.compile(r'.*\.gbk'))

for filename in filenames:
    rec_name = filename[:filename.find(".gbk")]
    print '.',

    # load data
    record = load_genbank(dir_in+"/"+filename)

    # scan annotations
    for feat in record.features:
        if feat.type == feat_type:
            try:
                if feat_name in feat.qualifiers[feat_tag]:
                    print '\nfound', feat_name, 'in', rec_name
                    # extract sequence
                    new_rec = feat.extract(record)
                    new_rec.id = rec_name+'_'+feat_name
                    new_rec.description = "Extracted from "+new_rec.description
                    records.append(new_rec)

            except KeyError:
                pass

print ''
write_fasta(main_out, records)
