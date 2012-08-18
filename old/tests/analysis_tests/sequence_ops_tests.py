__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from analysis import sequence_ops, seqfile_ops

class test_seqfeat_ops(TestCase):

    def setUp(self):
        # create temp directory
        self.temp_dir = "tests/temp_data/"
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
        self.filename = self.temp_dir+"temp_in.gbk"
        self.count = seqfile_ops.write_genbank(self.filename,
                                                    self.record)

    def tearDown(self):
        try: os.remove(self.filename)
        except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_feat_collect_all_cds(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feat_collect_all_genes(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('genes'), 'tags': {}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feat_collect_cds_by_locus_tag(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {'locus_tag': ['locustag 1',
                                                              'locustag 4']}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)
        
    def test_feat_collect_gene_by_locus_tag(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('gene'), 'tags': {'locus_tag': ['locustag 1',
                                                               'locustag 4']}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)

    def test_feat_collect_cds_by_product(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS'), 'tags': {'locus_tag': ['locustag 1',
                                                              'locustag 4']}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)

    def test_feat_collect_gene_by_product(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('gene'), 'tags': {'locus_tag': ['locustag 1',
                                                               'locustag 4']}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 1)

    def test_feat_collect_mixed(self):
        self.assertIs(self.count, 1)
        feat_mode = {'types': ('CDS', 'gene'), 'tags': {'locus_tag':
                                                            ('locustag 3'),
                                                        'product':
                                                            (('product 2'))}}
        collected = sequence_ops.feat_collect(self.filename, feat_mode)
        self.assertIs(len(collected), 2)

    def test_feature_coords(self):
        coords_list = sequence_ops.feature_coords(self.features)
        print coords_list
        self.assertIs(len(coords_list), len(self.features))
        index = 0
        for (coord1, coord2) in coords_list:
            self.assertEquals(str(self.features[index].location.start),
                          str(coord1))
            index +=1

class test_seqrecord_ops(TestCase):

    def setUp(self):
        self.sequence = Seq('AATTTAATGGCGCAGGCTAAGGCTCGTTTTTGGCGCT')
        self.record = SeqRecord(self.sequence, id='temp_seq')
        self.start = 5
        self.stop = 10
        self.coords_list = [(2,12),(15,30)]

    def tearDown(self):
        pass

    def test_get_region_seq(self):
        segment_seq = sequence_ops.get_region_seq(self.record, self.start, self.stop)
        self.assertEqual(str(segment_seq), 'AATGG')

    def test_get_seq_subset_by_coords(self):
        subset = sequence_ops.get_seq_subset_by_coords(self.record, self.coords_list)
        self.assertEqual(len(subset), len(self.coords_list))
        self.assertEqual(str(subset[0].seq), 'TTTAATGGCG')
        self.assertEqual(str(subset[1].seq), 'GCTAAGGCTCGTTTT')
        self.assertEqual(subset[0].id, 'temp_seq_2-12')
        self.assertEqual(subset[1].id, 'temp_seq_15-30')

class test_coord_chop(TestCase):

    def setUp(self):
        self.length = 7676
        self.size = 555
        self.count = 12
        self.coords_list = [(42,8604),(1534,30643)]

    def tearDown(self):
        pass

    def test_chop_exactsize(self):
        pair_list = sequence_ops.coord_chop(self.length, self.size, 'exact_size')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertEquals(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)
        self.assertGreaterEqual(first_pair_len, last_pair_len)

    def test_chop_maxsize_bisect(self):
        pair_list = sequence_ops.coord_chop(self.length, self.size, 'maxsize_bisect')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(abs(first_pair_len-last_pair_len), 1)

    def test_chop_maxsize_divisor(self):
        pair_list = sequence_ops.coord_chop(self.length, self.size, 'maxsize_divisor')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)

    def test_chop_count_divisor(self):
        pair_list = sequence_ops.coord_chop(self.length, self.count, 'count_divisor')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertLessEqual(abs(len(pair_list)-self.count), 1)
        self.assertGreaterEqual(first_pair_len, last_pair_len)

    def test_chop_none(self):
        pair_list = sequence_ops.coord_chop(self.length, None, None)
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        self.assertIs(len(pair_list), 1)
        self.assertEqual(first_pair_len, self.length)

    def test_recursive_bisector(self):
        pair_list = sequence_ops.recursive_bisector(self.coords_list, self.size)
        print pair_list
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)

        