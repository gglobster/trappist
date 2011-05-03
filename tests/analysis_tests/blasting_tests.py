__author__ = 'GG'

from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from analysis import sequence_file_ops, blasting

class test_basic_blasting(TestCase):

    def setUp(self): # move part of this to a fixtures file
        seq_1 = Seq('AATTTAATGGCGCAGGCT')
        seq_2 = Seq('GGATTATACCAAAGGCTTAAACTAT')
        seq_3 = Seq('GGAATATACCTTAGGCTTAAACTAT')
        record_1 = SeqRecord(seq_1, id='temp_1')
        record_2 = SeqRecord(seq_2, id='temp_2')
        record_3 = SeqRecord(seq_3, id='temp_3')
        seq_records_query = [record_1]
        seq_records_db = [record_1, record_2, record_3]
        self.query_file = 'temp_in.fas'
        sequence_file_ops.write_fasta(self.query_file, seq_records_query)
        db_file = 'temp_db.fas'
        sequence_file_ops.write_fasta(db_file, seq_records_db)
        self.db_name = 'genome'
        self.database = blasting.make_blastDB('', self.db_name, db_file,
                                              'nucl')
        self.out_file = 'temp_out.blast'
        self.prefs = {'evalue': 0.1, 'outfmt_pref': 6, 'score': 10,
                      'length': 10}

    def tearDown(self):
        pass  

    def test_local_blastn(self):
        self.status = blasting.local_blastn(self.query_file, self.out_file,
                                            self.db_name, self.prefs)
        self.assertIs(self.status['output'])
        self.assertIsNone(self.status['error'])
        self.matches = blasting.parse_blast_out6(self.out_file, self.prefs)
        
        