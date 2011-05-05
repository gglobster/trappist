__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio import AlignIO
from Bio.Seq import Seq, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from analysis import seq_features, sequence_file_ops

class test_seq_features(TestCase):

    def setUp(self):
        # create temp directory
        self.temp_dir = 'tests/temp_data/'
        os.mkdir(self.temp_dir)
        # create a sequence record with genbank features
        self.seq = Seq('AATTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGCGGGGCGGGG',
                       generic_dna)
        self.feat_cds1 = SeqFeature(FeatureLocation(5,10), strand=1,
                                    type='CDS', qualifiers={'locus_tag':
                                                                'locustag 1',
                                                            'product':
                                                                'product 1'})
        self.feat_cds2 = SeqFeature(FeatureLocation(22,40), strand=-1,
                                    type='CDS', qualifiers={'locus_tag':
                                                                'locustag 2',
                                                            'product':
                                                                'product 2'})
        self.feat_gene1 = SeqFeature(FeatureLocation(5,10), strand=-1,
                                     type='gene', qualifiers={'locus_tag':
                                                                'locustag 3',
                                                            'product':
                                                                'product 3'})
        self.feat_gene2 = SeqFeature(FeatureLocation(22,40), strand=-1,
                                     type='gene', qualifiers={'locus_tag':
                                                                'locustag 4',
                                                            'product':
                                                                'product 4'})
        self.features = [self.feat_cds1, self.feat_cds2, self.feat_gene1,
                         self.feat_gene2]
        self.record = SeqRecord(self.seq, id='temp', name='temp_name',
                                description='temporary genbank record',
                                features=self.features)
        self.filename = self.temp_dir+'temp_in.gbk'
        self.count = sequence_file_ops.write_genbank(self.filename,
                                                    self.record)

    def tearDown(self):
        try: os.remove(self.filename)
        except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_feat_collect_all_cds(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feat_collect_all_genes(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('genes'), 'tags': {}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feat_collect_cds_by_locus_tag(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {'locus_tag': ['locustag 1',
                                                              'locustag 4']}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)
        
    def test_feat_collect_gene_by_locus_tag(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('gene'), 'tags': {'locus_tag': ['locustag 1',
                                                               'locustag 4']}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)

    def test_feat_collect_cds_by_product(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {'locus_tag': ['locustag 1',
                                                              'locustag 4']}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)

    def test_feat_collect_gene_by_product(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('gene'), 'tags': {'locus_tag': ['locustag 1',
                                                               'locustag 4']}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)


    def test_feat_collect_mixed(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS', 'gene'), 'tags': {'locus_tag':
                                                            ('locustag 3'),
                                                        'product':
                                                            (('product 2'))}}
        collected = seq_features.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feature_coords(self):
        coords_list = seq_features.feature_coords(self.features)
        print coords_list
        self.assertIs(len(coords_list), len(self.features))
        index = 0
        for (coord1, coord2) in coords_list:
            self.assertEquals(str(self.features[index].location.start),
                          str(coord1))
            index +=1

    