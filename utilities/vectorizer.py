# script to represent sequences as vectors (based on gene annots)

from sys import argv
from libs.common import load_genbank, write_fasta, make_blastDB, \
    local_tblastn_2file, read_array, blast_dtypes, ensure_dir
from libs.tetris import segment_finder
from Bio.SeqRecord import SeqRecord
import cPickle as pickle

from sets.NheA_ctxt_set import test as genomes

data_dir = 'data/'+argv[1]+'/'
seq_dir = data_dir+argv[2]+'/'
out_dir = data_dir+argv[3]+'/'
feat_type = argv[4]
threshold = int(argv[5])
min_com = int(argv[6]) # min number of non-core feats in common within groups

ensure_dir([out_dir])

db_file = out_dir+'ref_DB.fas'
db_path = out_dir+'refs'

core_genome_pickle = out_dir+'core_genome.pik'
cluster_set_file = out_dir+'clusters.py'

new_DB = True
init_DB = False

symbolDB = {}
vectorDB = {}
segmentDB = {}

# vectorDB is a dict that contains genome-keyed dicts,
# where each key-value pair is a version of the vectorized genomes:
# 'symbols' is the basic sequence of symbols that makes up the vector;
# 'red_vector' is the primary reduced form obtained by grouping syntenous
# core symbols ;

core_keys = []

sym_cnt = 0

# check for pickled info
try:
    pickle_guts = pickle.load(open(core_genome_pickle, 'rb'))
    core_keys = pickle_guts['core']
    segmentDB = pickle_guts['core_segs']
    vectorDB = pickle_guts['vectors']

except IOError: redo_flag = True
except KeyError: redo_flag = True
else: redo_flag = False

if redo_flag:

    ### PART 1 : Vectorize ###

    for genome in genomes:

        g_vector = []

        print genome['name'],

        while True:

            try:
                assert genome['input'] == 'gbk'
            except ValueError:
                print "bad format (skipping)"
                break

            # load genome file to extract features (to proteins in mfas file)
            record = load_genbank(seq_dir+genome['file'])
            select = [feat for feat in record.features
                      if feat.type == feat_type]
            feat_cnt = 0

            # cycle through selected features
            for feat in select:

                feat_cnt +=1
                rec = feat.extract(record)
                rec.description = genome['name']+'_'+feat_type+'_'+str(feat_cnt)

                # initialize or update blast DB
                if init_DB:
                    ref_records = [value[0] for value in symbolDB.values()]
                    write_fasta(db_file, ref_records)
                    try:
                        make_blastDB(db_path, db_file, 'nucl')
                    except Exception:
                        print "failed to make blast DB"
                        exit()
                    init_DB = False

                # first go: add all features as new symbols
                if new_DB:
                    sym_cnt +=1
                    symbol = 'N'+str(sym_cnt)
                    rec.id = symbol
                    symbolDB[symbol] = [rec]
                    g_vector.append(symbol)

                else:
                    # tblastn against the reference DB
                    infile = data_dir+'temp.fas'
                    outfile = data_dir+'temp.txt'
                    write_fasta(infile, SeqRecord(rec.seq.translate(), id='temp'))
                    prefs = {'evalue': 0.001, 'outfmt_pref': 6}
                    try:
                        local_tblastn_2file(infile, db_path, outfile, prefs)
                    except Exception:
                        print "failed to blast"
                        exit()

                    # parse output -- take only first hit
                    try:
                        hit = read_array(outfile, blast_dtypes)[0]
                    except IndexError:
                        print '-',
                        # add record to DB as new symbol
                        sym_cnt +=1
                        symbol = 'N'+str(sym_cnt)
                        rec.id = symbol
                        symbolDB[symbol] = [rec]
                        g_vector.append(symbol)
                        init_DB = True
                    else:
                        if hit[2] > threshold:
                            print '+',
                            # add record to DB as secondary match
                            symbol = hit[1]
                            rec.id = symbol+'_'+str(len(symbolDB[symbol])+1)
                            symbolDB[symbol].append(rec)
                            # add symbol to genome vector
                            g_vector.append(symbol)
                            init_DB = False
                        else:
                            print '<',
                            # add record to DB as new symbol
                            sym_cnt +=1
                            symbol = 'N'+str(sym_cnt)
                            rec.id = symbol
                            symbolDB[symbol] = [rec]
                            g_vector.append(symbol)
                            init_DB = True

            break

        vectorDB[genome['name']] = {'symbols': g_vector}

        if new_DB:
            init_DB = True
            new_DB = False

        print 'OK'

    # find core genome
    core_genome = [value for value in symbolDB.values()
                   if len(value) == len(genomes)]
    core_keys = [key for key in symbolDB
                 if len(symbolDB[key]) == len(genomes)]

    # find core segments (only need one vector as ref to get universal segs)
    seg_group_dict = {}
    init_seg_group = True
    reference = [vectorDB[v_key]['symbols'] for v_key in vectorDB][0]
    segmentDB = segment_finder(reference, vectorDB, 'symbols', 'SG')
    core_segment_symbols = []
    solon_cnt = 0
    for segment in segmentDB.values():
        for symbol in segment:
            core_segment_symbols.append(symbol)
    for core_key in core_keys:
        if core_key not in core_segment_symbols:
            segmentDB['SL'+str(solon_cnt+1)] = [core_key]

    # refresh core segment symbols
    core_segment_symbols = []
    for segment in segmentDB.values():
        for symbol in segment:
            core_segment_symbols.append(symbol)

    # generate reduced vectors (with core segments)
    for v_key in vectorDB:
        red_vector = []
        for symbol in vectorDB[v_key]['symbols']:
            if symbol not in core_segment_symbols:
                red_vector.append(symbol)
            else:
                try:
                    red_vector.append([seg_group for seg_group in segmentDB if
                                       symbol in segmentDB[seg_group] and
                                       seg_group not in red_vector][0])
                except IndexError:
                    pass
        vectorDB[v_key]['red_vector'] = red_vector

    # output protein sequences by feature
    for nt_set in core_genome:
        aa_set = [SeqRecord(nt_rec.seq.translate(), id=nt_rec.id+'_aa',
                            description=nt_rec.description)
                  for nt_rec in nt_set]
        write_fasta(out_dir+nt_set[0].id+'_aa.fas', aa_set)

    # concatenate aa seqs for later phylo analysis
    counter = 0
    concat_recs = []
    while counter < len(genomes):
        nt_recs = [sym_set[counter] for sym_set in core_genome]
        aa_recs = [SeqRecord(nt_rec.seq.translate(), id=nt_rec.id+'_aa',
                             description=nt_rec.description)
                   for nt_rec in nt_recs]

        aa_concat = SeqRecord('', '')
        for aa_rec in aa_recs:
            aa_concat += aa_rec

        aa_concat.id = genomes[counter]['name']
        aa_concat.description = ':'.join(core_keys)

        concat_recs.append(aa_concat)
        counter +=1

    # write concatenated sequences to file
    write_fasta(out_dir+'concat_aa.fas', concat_recs)

    # save core and vectors in a pickle
    pickle_guts = {'core': core_keys,
                   'core_segs': segmentDB,
                   'vectors': vectorDB}
    pickle.dump(pickle_guts, open(core_genome_pickle, 'wb'))

else:
    print "Loaded pickled core genome and vectors"

print "Core genome:", len(core_keys), "features in",\
                      len(segmentDB), "segments"

print segmentDB

### PART 2 : Cluster ###

# group vectors by min common (besides core genome)
clusters_dict = {}
init_cluster = True

for v_key in vectorDB:

    if init_cluster:
        clusters_dict['C'+str(len(clusters_dict)+1)] = {
            'commons': vectorDB[v_key]['red_vector'], 'genomes': [v_key]}
        init_cluster = False

    else:
        vector_assigned = False
        for cluster in clusters_dict:
            new_commons = []
            for item in vectorDB[v_key]['red_vector']:
                if item in clusters_dict[cluster]['commons']:
                    new_commons.append(item)
            if len(new_commons) >= len(core_keys)+min_com:
                vector_assigned = True
                clusters_dict[cluster]['commons'] = new_commons
                clusters_dict[cluster]['genomes'].append(v_key)

        if not vector_assigned:
            clusters_dict['C'+str(len(clusters_dict)+1)] = {
                'commons': vectorDB[v_key]['red_vector'], 'genomes': [v_key]}

print "Found", len(clusters_dict), "clusters with at least", \
      len(core_keys)+min_com, "features in common"

# find conserved segments within clusters
for cluster in clusters_dict:
    reference = clusters_dict[cluster]['commons']
    c_vectorDB = {}
    for v_key in clusters_dict[cluster]['genomes']:
        c_vectorDB[v_key] = vectorDB[v_key]
    cluster_segmentDB = segment_finder(reference, c_vectorDB, 'red_vector',
                                       'CG')
    print cluster_segmentDB

# output as python-loadable genome sets
set_lines = []
name_check = {}
for cluster in clusters_dict:
    set_lines.append("".join([cluster, " = ["]))
    order = 0
    for g_name in clusters_dict[cluster]['genomes']:
        if g_name in name_check:
            name_check[g_name].append(cluster)
        else:
            name_check[g_name] = [cluster]
        genome = [genome for genome in genomes if genome['name'] == g_name][0]
        order +=1
        line = "".join(["\t{'name': '", g_name,
                        "', 'file': '", genome['file'],
                        "', 'input': '", genome['input'],
                        "', 'order': ", str(order),
                        ", 'nudge': 0, 'offset': 0,'ignore': (0, 0)},"])
        set_lines.append(line)
    set_lines.append("]\n")

open(cluster_set_file, 'w').write("\n".join(set_lines))

# check for overlaps (genomes assigned several clusters)
print "Overlap alerts:"
counter = 0
for g_name in name_check:
    if len(name_check[g_name]) > 1:
        print "\t", g_name, name_check[g_name]
        counter +=1
if not counter:
    print "\tnone"
