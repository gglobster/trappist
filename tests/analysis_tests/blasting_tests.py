__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from analysis import sequence_file_ops, blasting

class test_blasting(TestCase):

    def setUp(self): # move part of this to a fixtures file
        self.temp_dir = 'tests/temp_data/'
        os.mkdir(self.temp_dir)
        seq_1 = Seq('AATTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGCGGCGGGCGCGGCG')
        seq_2 = Seq('GGATTATACCAAAGGCTTAAACTATAGGCTAGGAGAGATAGACG')
        seq_3 = Seq('GGAATATACCTTAGGCTTAAACTATAGGCTAGGAGAGGCTCG')
        record_1 = SeqRecord(seq_1, id='temp_1')
        record_2 = SeqRecord(seq_2, id='temp_2')
        record_3 = SeqRecord(seq_3, id='temp_3')
        seq_records_query = [record_1]
        seq_records_db = [record_1, record_2, record_3]
        self.query_file = self.temp_dir+'temp_in.fas'
        sequence_file_ops.write_fasta(self.query_file, seq_records_query)
        self.db_file = self.temp_dir+'temp_db.fas'
        sequence_file_ops.write_fasta(self.db_file, seq_records_db)
        self.db_name = 'genome'
        self.database = blasting.make_blastDB(self.temp_dir, self.db_name,
                                              self.db_file, 'nucl')
        self.out_file = self.temp_dir+'temp_out.blast'
        self.prefs = {'evalue': 0.1, 'outfmt_pref': 6, 'score': 10,
                      'length': 10}

    def tearDown(self):
        temp_files = [self.query_file, self.db_file, self.out_file]
        for file in temp_files:
            os.remove(file)
        blast_db_ext_set = ['.nin', '.nog', '.nsd', '.nsi', '.nsq', '.nhr']
        for blast_db_ext in blast_db_ext_set:
            os.remove(self.temp_dir+self.db_name+blast_db_ext)
        os.rmdir(self.temp_dir)

    def test_local_blastn(self):
        self.status = blasting.local_blastn(self.query_file, self.out_file,
                                            self.temp_dir+self.db_name,
                                            self.prefs)
        self.assertEquals(self.status['output'], '')
        self.assertIsNone(self.status['error'])
        self.matches = blasting.parse_blast_out6(self.out_file, self.prefs)
        self.assertEquals(len(self.matches), 1)
        
        
        