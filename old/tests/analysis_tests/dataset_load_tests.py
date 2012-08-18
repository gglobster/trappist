__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from analysis import seqfile_ops, dataset_load

class test_seq_subset_load(TestCase):

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
        # create some fasta sequence records
        self.seq_1 = Seq('AATTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGCGGCGGGG')
        self.seq_2 = Seq('GGATTATACCAAAGGCTTAAACTATAGGCTAGGAGAGATAGACG')
        self.seq_3 = Seq('GGAATATACCTTAGGCTTAAACTATAGGCTAGGAGAGGCTCG')
        self.seq_4 = Seq('GGGGATTACAGCCATAGTAACCAGATATTAAGACG')
        self.seq_5 = Seq('GGAACCGCTGATACATGATTATAGATCTATAGGGTCTAAAACATCG')
        self.record_1 = SeqRecord(self.seq_1, id='temp_1')
        self.record_2 = SeqRecord(self.seq_2, id='temp_2')
        self.record_3 = SeqRecord(self.seq_3, id='temp_3')
        self.record_4 = SeqRecord(self.seq_4, id='temp_4')
        self.record_5 = SeqRecord(self.seq_5, id='temp_5')
        # create record sets
        self.single_record = self.record_1
        self.multi_records = [self.record_1, self.record_2, self.record_3,
                              self.record_4, self.record_5]
        # create content for the coords file
        self.header_contents = ['pairs', 'start', 'stop']
        self.line_1_contents = ['pair_1', '0', '10']
        self.line_2_contents = ['pair_2', '12', '25']
        self.line_3_contents = ['pair_3', '30', '48']
        self.raw_contents = [self.header_contents, self.line_1_contents,
                             self.line_2_contents, self.line_3_contents]
        # transform into a single string
        self.str_contents = ""
        for line_set in self.raw_contents:
            joined = "\t".join(line_set)
            self.str_contents += joined+"\n"
        # define file names
        self.gbk_filename = self.temp_dir+"temp_in.gbk"
        self.single_q_file = self.temp_dir+"temp_single_in.fas"
        self.multi_q_file = self.temp_dir+"temp_multi_in.fas"
        self.coords_file = self.temp_dir+"temp_coords.txt"

    def tearDown(self):
        temp_files = [self.gbk_filename, self.single_q_file,
                      self.multi_q_file, self.coords_file]
        for file in temp_files:
            try: os.remove(file)
            except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_seq_subset_load_from_multifasta(self):
        seqfile_ops.write_fasta(self.multi_q_file, self.multi_records)
        subset_mode = 'flatfile'
        subset_args = None
        subset, subset_file = dataset_load.seq_subset_load(self.multi_q_file,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 5)
        index = 0
        for record in subset:
            self.assertEqual(subset[index].id, self.multi_records[index].id)
            index += 1
        self.assertIs(subset_file, self.multi_q_file)

    def test_seq_subset_load_from_coords(self):
        seqfile_ops.write_fasta(self.single_q_file, self.single_record)
        temp_file = open(self.coords_file, 'w')
        temp_file.write(self.str_contents)
        temp_file.close()
        subset_mode = 'coordinates'
        subset_args = {'file': self.coords_file, 'header': 1,
                       'columns': (1, 2)}
        subset, subset_file = dataset_load.seq_subset_load(self.single_q_file,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 3)
        self.assertEqual(subset[0].id, 'temp_1_0-10')
        self.assertEqual(str(subset[2].seq), 'TTTGGCGCTCGCGGCGGG')

    def test_seq_subset_load_from_all_cds(self):
        seqfile_ops.write_genbank(self.gbk_filename, self.record)
        subset_mode = 'features'
        subset_args = {'types': ('CDS'), 'tags': {}}
        subset, subset_file = dataset_load.seq_subset_load(self.gbk_filename,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 2)

    def test_seq_subset_load_from_all_genes(self):
        seqfile_ops.write_genbank(self.gbk_filename, self.record)
        subset_mode = 'features'
        subset_args = {'types': ('genes'), 'tags': {}}
        subset, subset_file = dataset_load.seq_subset_load(self.gbk_filename,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 2)

    def test_seq_subset_load_from_gene_by_product(self):
        seqfile_ops.write_genbank(self.gbk_filename, self.record)
        subset_mode = 'features'
        subset_args = {'types': ('gene'), 'tags': {'locus_tag': ['locustag 1',
                                                               'locustag 4']}}
        subset, subset_file = dataset_load.seq_subset_load(self.gbk_filename,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 1)
        self.assertEqual(subset[0].id, 'temp_22-40')

    def test_seq_subset_load_from_cds_by_locus_tag(self):
        seqfile_ops.write_genbank(self.gbk_filename, self.record)
        subset_mode = 'features'
        subset_args = {'types': ('CDS'), 'tags': {'locus_tag': ['locustag 1',
                                                              'locustag 4']}}
        subset, subset_file = dataset_load.seq_subset_load(self.gbk_filename,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 1)
        self.assertEqual(subset[0].id, 'temp_5-10')

    def test_seq_subset_load_from_mixed_features(self):
        seqfile_ops.write_genbank(self.gbk_filename, self.record)
        subset_mode = 'features'
        subset_args = {'types': ('CDS', 'gene'), 'tags': {'locus_tag':
                                                            ('locustag 3'),
                                                        'product':
                                                            (('product 2'))}}
        subset, subset_file = dataset_load.seq_subset_load(self.gbk_filename,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 2)
 
    def test_seq_subset_load_from_chop_by_size(self):
        seqfile_ops.write_fasta(self.single_q_file, self.single_record)
        subset_mode = 'size'
        subset_args = {'size': 5, 'chop_mode': 'exact_size'}
        subset, subset_file = dataset_load.seq_subset_load(self.single_q_file,
                                                           subset_mode,
                                                           subset_args)
        self.assertIs(len(subset), 10)
        self.assertEqual(subset[0].id, 'temp_1_0-5')

    # TODO: could add tests for other chop modes that I'm too lazy to do now


class test_genome_sets_load(TestCase):

    def setUp(self):
        # create temp directory
        self.temp_dir = "tests/temp_data/"
        os.mkdir(self.temp_dir)
        # assign database directory (but don't create it yet)
        self.db_path = "tests/temp_data/temp_db/"
        # define file names
        self.gen_filename_1 = "temp_gen1.fas"
        self.gen_filename_2 = "temp_gen2.fas"
        self.gen_filename_3 = "temp_gen3.fas"
        self.genomes_list = self.temp_dir+"temp_genomes.txt"
        # create content for the genomes file
        self.header_contents = ['genome_name', 'file']
        self.line_1_contents = ['genome_1', self.gen_filename_1]
        self.line_2_contents = ['genome_2', self.gen_filename_2]
        self.line_3_contents = ['genome_3', self.gen_filename_3]
        self.raw_contents = [self.header_contents, self.line_1_contents,
                             self.line_2_contents, self.line_3_contents]
        # transform into a single string
        self.str_contents = ""
        for line_set in self.raw_contents:
            joined = "\t".join(line_set)
            self.str_contents += joined+"\n"
        # create some sequence records
        self.seq_1 = Seq('AATTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGCGGCGGGG')
        self.seq_2 = Seq('GGATTATACCAAAGGCTTAAACTATAGGCTAGGAGAGATAGACG')
        self.seq_3 = Seq('GGAATATACCTTAGGCTTAAACTATAGGCTAGGAGAGGCTCG')
        self.seq_4 = Seq('GGGGATTACAGCCATAGTAACCAGATATTAaGACG')
        self.seq_5 = Seq('GGAACCGCTGATACATGATTATAGATCTATAGGGTCTAAAACATCG')
        self.seq_6 = Seq('AGGTCATGTACGATGCAGAATTTGTCGTACGATGTTAGTACGATGGTA')
        self.seq_7 = Seq('TTTTTCGCGCGCTTAGACCCAAAATATATTGTCGCTATAGGTCCCTCT')
        self.seq_8 = Seq('ACCGTGTGGCATTTATATTACACCACACACAGATTGGGTGTGCCAATCAG')
        self.seq_9 = Seq('ACCGTACGTACCATATTATTATATAGGATAGATATTTAGAGGATTTAGAT')
        self.record_1 = SeqRecord(self.seq_1, id='temp_1')
        self.record_2 = SeqRecord(self.seq_2, id='temp_2')
        self.record_3 = SeqRecord(self.seq_3, id='temp_3')
        self.record_4 = SeqRecord(self.seq_4, id='temp_4')
        self.record_5 = SeqRecord(self.seq_5, id='temp_5')
        self.record_6 = SeqRecord(self.seq_5, id='temp_6')
        self.record_7 = SeqRecord(self.seq_5, id='temp_7')
        self.record_8 = SeqRecord(self.seq_5, id='temp_8')
        self.record_9 = SeqRecord(self.seq_5, id='temp_9')
        # create record sets
        self.gen1_records = [self.record_1, self.record_2, self.record_3]
        self.gen2_records = [self.record_4, self.record_5, self.record_6]
        self.gen3_records = [self.record_7, self.record_8, self.record_9]
        # create all primary files
        seqfile_ops.write_fasta(self.temp_dir+self.gen_filename_1,
                                self.gen1_records)
        seqfile_ops.write_fasta(self.temp_dir+self.gen_filename_2,
                                self.gen2_records)
        seqfile_ops.write_fasta(self.temp_dir+self.gen_filename_3,
                                self.gen3_records) 
        seq_file = open(self.genomes_list, 'w')
        seq_file.write(self.str_contents)
        seq_file.close()

    def tearDown(self):
        temp_files = [self.temp_dir+self.gen_filename_1,
                      self.temp_dir+self.gen_filename_2,
                      self.temp_dir+self.gen_filename_3,
                      self.genomes_list]
        for file in temp_files:
            try: os.remove(file)
            except Exception as message: print message
        blast_db_ext_set = [".nin", ".nog", ".nsd", ".nsi", ".nsq", ".nhr"]
        genome_list = ['genome_1', 'genome_2', 'genome_3']
        for genome in genome_list:
            for blast_db_ext in blast_db_ext_set:
                try: os.remove(self.db_path+genome+blast_db_ext)
                except Exception as message: print message
        try: os.rmdir(self.db_path)
        except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_genome_sets_load(self):
        input_prefs = {'header': 1, 'columns': (0, 1)}
        genome_sets = dataset_load.genome_sets_load(self.temp_dir,
                                                    self.genomes_list,
                                                    input_prefs,
                                                    self.db_path)
        self.assertIs(len(genome_sets), 3)
        self.assertEqual(genome_sets[0].name, 'genome_1')
        self.assertEqual(genome_sets[1].name, 'genome_2')
        self.assertEqual(genome_sets[2].name, 'genome_3')
        