__author__ = 'GG'

# TODO: file format tests could be more complete

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from analysis import seqfile_ops

class test_seqfile_ops(TestCase):

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
        self.feat_gene1 = SeqFeature(FeatureLocation(5,10), strand=-1,
                                     type='gene', qualifiers={'locus_tag':
                                                                'locustag 3',
                                                            'product':
                                                                'product 3'})
        self.features = [self.feat_cds1, self.feat_gene1]
        self.record = SeqRecord(self.seq, id='temp', name='temp_name',
                                description='temporary genbank record',
                                features=self.features)
        # create some extra sequence records for the multifasta file
        self.seq_1 = Seq('GGTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGAAGCGG')
        self.seq_2 = Seq('GGTTAATGCCGCAGGCTAGGAGAGAGATGTGGCGCTCGCGGCGGCG')
        self.record_1 = SeqRecord(self.seq_1, id='temp_1')
        self.record_2 = SeqRecord(self.seq_2, id='temp_2')
        # create a set of records for the multifasta file
        self.three_records = [self.record, self.record_1, self.record_2]
        # define file names
        self.gbk_filename = self.temp_dir+"temp_in.gbk"
        self.fas_filename = self.temp_dir+"temp_in.fas"
        self.seq_filename = self.temp_dir+"temp_in.seq"
        self.other_filename = self.temp_dir+"temp_in.txt"

    def tearDown(self):
        temp_files = [self.gbk_filename, self.fas_filename,
                      self.seq_filename, self.other_filename]
        for file in temp_files:
            try: os.remove(file)
            except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_seqfile_format_fas(self):
        count = seqfile_ops.write_fasta(self.fas_filename, self.record)
        self.assertIs(count, 1)
        format, name = seqfile_ops.seqfile_format(self.fas_filename)
        self.assertEqual(format, 'fasta')

    def test_seqfile_format_gbk(self):
        count = seqfile_ops.write_genbank(self.gbk_filename,
                                                self.record)
        self.assertIs(count, 1)
        format, name = seqfile_ops.seqfile_format(self.gbk_filename)
        self.assertEqual(format, 'genbank')

    def test_seqfile_format_seq(self):
        seq_file = open(self.seq_filename, 'w')
        seq_file.write(str(self.seq))
        seq_file.close()
        format, name = seqfile_ops.seqfile_format(self.seq_filename)
        self.assertEqual(format, 'seq')

    def test_seqfile_format_other(self):
        seq_file = open(self.other_filename, 'w')
        seq_file.write(str(self.seq))
        seq_file.close()
        format, name = seqfile_ops.seqfile_format(self.other_filename)
        self.assertEqual(format, 'unidentified')

    def test_surefmt_load_fas2fas(self):
        count = seqfile_ops.write_fasta(self.fas_filename, self.record)
        self.assertIs(count, 1)
        fas_record = seqfile_ops.surefmt_load(self.fas_filename,
                                                    'fasta', generic_dna)
        self.assertEqual(fas_record.id, self.record.id)

    def test_surefmt_load_gbk2fas(self):
        count = seqfile_ops.write_genbank(self.gbk_filename,
                                                self.record)
        self.assertIs(count, 1)
        fas_record = seqfile_ops.surefmt_load(self.gbk_filename,
                                                    'fasta', generic_dna)
        self.assertEqual(fas_record.id, self.record.id)

    def test_surefmt_load_gbk2gbk(self):
        count = seqfile_ops.write_genbank(self.gbk_filename,
                                                self.record)
        self.assertIs(count, 1)
        gbk_record = seqfile_ops.surefmt_load(self.gbk_filename,
                                                    'genbank', generic_dna)
        self.assertEqual(gbk_record.id, self.record.id)
        # check features
        for index in range (0,1):
            self.assertEqual(gbk_record.features[index].type,
                             self.record.features[index].type)

    def test_surefmt_load_fas2gbk(self):
        count = seqfile_ops.write_fasta(self.fas_filename, self.record)
        self.assertIs(count, 1)
        gbk_record = seqfile_ops.surefmt_load(self.fas_filename,
                                                    'genbank', generic_dna)
        self.assertEqual(gbk_record.id, self.record.id)

    def test_write_and_load_single_fasta(self):
        count = seqfile_ops.write_fasta(self.fas_filename, self.record)
        self.assertIs(count, 1)
        fas_record = seqfile_ops.load_fasta(self.fas_filename)
        self.assertEqual(fas_record.id, self.record.id)

    def test_write_and_load_multifasta(self):
        count = seqfile_ops.write_fasta(self.fas_filename,
                                              self.three_records) 
        self.assertIs(count, 3)
        fas_records = seqfile_ops.load_multifasta(self.fas_filename)
        for index in range (0,2):
            self.assertEqual(fas_records[index].id,
                             self.three_records[index].id)

    def test_write_and_load_single_genbank(self):
        count = seqfile_ops.write_genbank(self.gbk_filename,
                                                self.record)
        self.assertIs(count, 1)
        gbk_record = seqfile_ops.load_genbank(self.gbk_filename)
        self.assertEqual(gbk_record.id, self.record.id)
        # check features
        for index in range (0,1):
            self.assertEqual(gbk_record.features[index].type,
                             self.record.features[index].type)

    # TODO: write test_load_multiple_genbank

    def test_load_agnostic_fas(self):
        count = seqfile_ops.write_fasta(self.fas_filename, self.record)
        self.assertIs(count, 1)
        fas_record, type = seqfile_ops.load_agnostic(self.fas_filename)
        self.assertEqual(fas_record.id, self.record.id)
        self.assertEqual(type, 'fasta')

    def test_load_agnostic_gbk(self):
        count = seqfile_ops.write_genbank(self.gbk_filename,
                                                self.record)
        self.assertIs(count, 1)
        gbk_record, type = seqfile_ops.load_agnostic(self.gbk_filename)
        self.assertEqual(gbk_record.id, self.record.id)
        self.assertEqual(type, 'genbank')
        # check features
        for index in range (0,1):
            self.assertEqual(gbk_record.features[index].type,
                             self.record.features[index].type)

    # TODO: write test_compile_multifasta_from_folders(self):


