from array_tetris import coord_chop
from Bio.SeqFeature import SeqFeature, FeatureLocation
from writers import write_fasta
from loaders import load_genbank

class Reference(object):
    """Persistent reference object."""

    def __init__(self, name, file, input, seg_mode, capture, fas_out, gbk_out,
                 seg_out, logfile):
        self.name = name
        self.file = file
        self.input = input
        self.seg_mode = seg_mode
        self.segs = []
        self.capture = capture
        self.fas = fas_out
        self.gbk = gbk_out
        self.segs_dir = seg_out
        self.logfile = logfile

    def get_segs_from_list(self, list):
        self.segs = list

    def get_segs_from_chop(self, length, size):
        pair_list = coord_chop(length, size, 'exact_size')
        counter = 0
        for a,b in pair_list:
            counter +=1
            seg = {'coords': (a, b), 'strand': 1, 'name': str(counter),
                   'note': str(a)+'_'+str(b)}
            self.segs.append(seg)

    def get_segs_from_feats(self, feat_type):
        feats = [feat for feat in load_genbank(self.gbk).features
                 if feat.type == feat_type]
        counter = 0
        for feat in feats: # TODO: there must be a better way to do this !!!
            counter +=1
            a = int(str(feat.location.start))
            b = int(str(feat.location.end))
            feat_id = feat_type+'_'+str(counter)
            seg = {'coords': (a, b), 'strand': feat.strand, 'name': feat_id,
                   'note': str(a)+'_'+str(b)}
            self.segs.append(seg)

    def extract_segs_seqs(self, record, out_dir):
        count = 0
        for seg in self.segs:
            # unpack segment coords
            seg_start, seg_stop = seg['coords'][0], seg['coords'][1]
            # extract segment sequence
            segment = record[seg_start:seg_stop]
            if seg['strand'] < 0:
                segment = segment.reverse_complement()
            segment.id = self.name+"_"+seg['name']
            # write to individual file
            out_file = out_dir+self.name+"_"+seg['name']+".fas"
            write_fasta(out_file, segment)
            # record segment feature
            feat_loc = FeatureLocation(seg_start, seg_stop)
            feature = SeqFeature(location=feat_loc,
                                 type='ref_seg',
                                 qualifiers={'id': seg['name']})
            record.features.append(feature)
            count +=1
        return record

    def log(self, string):
        ref_log = open(self.logfile, 'a')
        ref_log.write(string)
        ref_log.close()
