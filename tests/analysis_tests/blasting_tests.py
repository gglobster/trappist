__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from analysis import seqfile_ops, blasting

class test_blast_calls(TestCase):

    def setUp(self): # could move part of this to a fixtures file for reuse
        # create temp directory
        self.temp_dir = "tests/temp_data/"
        os.mkdir(self.temp_dir)
        # create some sequence records
        self.seq_1 = Seq('AATTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGCGGCGGGG')
        self.seq_2 = Seq('GGATTATACCAAAGGCTTAAACTATAGGCTAGGAGAGATAGACG')
        self.seq_3 = Seq('GGAATATACCTTAGGCTTAAACTATAGGCTAGGAGAGGCTCG')
        self.seq_4 = Seq('GGGGATTACAGCCATAGTAACCAGATATTAaGACG')
        self.seq_5 = Seq('GGAACCGCTGATACATGATTATAGATCTATAGGGTCTAAAACATCG')
        self.record_1 = SeqRecord(self.seq_1, id='temp_1')
        self.record_2 = SeqRecord(self.seq_2, id='temp_2')
        self.record_3 = SeqRecord(self.seq_3, id='temp_3')
        self.record_4 = SeqRecord(self.seq_4, id='temp_4')
        self.record_5 = SeqRecord(self.seq_5, id='temp_5')
        # create record sets
        self.single_record = self.record_1
        self.multi_records = [self.record_2, self.record_3, self.record_4]
        self.db_records = [self.record_1, self.record_2, self.record_3,
                           self.record_4, self.record_5]
        # define file names
        self.single_q_file = self.temp_dir+"temp_single_in.fas"
        self.multi_q_file = self.temp_dir+"temp_multi_in.fas"
        self.db_file = self.temp_dir+"temp_db.fas"
        self.db_name = 'temp_database'
        self.single_out_file = self.temp_dir+"temp_single_out.blast"
        self.multi_out_file = self.temp_dir+"temp_multi_out.blast"
        # define blast parameters
        self.prefs = {'evalue': 0.1, 'outfmt_pref': 6, 'score': 10,
                      'length': 10}

    def tearDown(self):
        temp_files = [self.single_q_file, self.multi_q_file, self.db_file,
                      self.single_out_file, self.multi_out_file]
        for file in temp_files:
            try: os.remove(file)
            except Exception as message: print message
        blast_db_ext_set = [".nin", ".nog", ".nsd", ".nsi", ".nsq", ".nhr"]
        for blast_db_ext in blast_db_ext_set:
            try: os.remove(self.temp_dir+self.db_name+blast_db_ext)
            except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_local_blastn(self):
        # prepare query
        seqfile_ops.write_fasta(self.single_q_file, self.single_record)
        query_record = seqfile_ops.load_fasta(self.single_q_file)
        self.assertEqual(query_record.id,self.record_1.id)
        self.assertEqual(str(query_record.seq),str(self.record_1.seq))
        # prepare database
        seqfile_ops.write_fasta(self.db_file, self.db_records)
        records_list = seqfile_ops.load_multifasta(self.db_file)
        index = 0
        for record in records_list:
            self.assertEqual(record.id,self.db_records[index].id)
            self.assertEqual(str(record.seq),str(self.db_records[index].seq))
            index +=1
        # make database
        self.dbfile_path, db_report = blasting.make_blastDB(self.temp_dir,
                                                            self.db_name,
                                                            self.db_file,
                                                            'nucl')
        self.assertIs(db_report['status'], 0)
        self.assertEquals(db_report['message'], 'database exists')
        # run local blast with single query
        self.status = blasting.local_blastn(self.single_q_file,
                                            self.single_out_file,
                                            self.dbfile_path,
                                            self.prefs)
        self.assertEquals(self.status['output'], '')
        self.assertIsNone(self.status['error'])
        # parse blast output
        matches_single = blasting.parse_blast_out6(self.single_out_file,
                                                   self.prefs)
        self.assertIs(len(matches_single), 1)
        self.assertEqual(matches_single[0]['contig_id'],
                         self.single_record.id)
        self.assertEqual(matches_single[0]['details']['match_p100'], 100)
        
    def test_blast_record_set(self):
        # prepare database
        seqfile_ops.write_fasta(self.db_file, self.db_records)
        db_records_list = seqfile_ops.load_multifasta(self.db_file)
        index = 0
        for record in db_records_list:
            self.assertEqual(record.id,self.db_records[index].id)
            self.assertEqual(str(record.seq),str(self.db_records[index].seq))
            index +=1
        # make database
        self.dbfile_path, db_report = blasting.make_blastDB(self.temp_dir,
                                                            self.db_name,
                                                            self.db_file,
                                                            'nucl')
        self.assertIs(db_report['status'], 0)
        self.assertEquals(db_report['message'], 'database exists')
        # run local blast batch (with multiple queries)
        matches_multi = blasting.blast_record_set(self.dbfile_path,
                                                  self.multi_records,
                                                  self.prefs)
        self.assertIs(len(matches_multi), 3)
        index = 0
        for record in self.multi_records:
            self.assertEqual(matches_multi[record.id][0]['contig_id'],
                             self.multi_records[index].id)
            self.assertEqual(matches_multi[record
            .id][0]['details']['match_p100'], 100) 
            index +=1



        
        