__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio.Seq import Seq, SeqRecord
from analysis import other

class test_seqrecord_stuff(TestCase):

    def setUp(self):
        self.sequence = Seq('AATTTAATGGCGCAGGCTAAGGCTCGTTTTTGGCGCT')
        self.record = SeqRecord(self.sequence, id='temp_seq')
        self.start = 5
        self.stop = 10
        self.coords_list = [(2,12),(15,30)]
    
    def tearDown(self):
        pass

    def test_get_region_seq(self):
        segment_seq = other.get_region_seq(self.record, self.start, self.stop)
        self.assertEqual(str(segment_seq), 'AATGG')

    def test_get_seq_subset_by_coords(self):
        subset = other.get_seq_subset_by_coords(self.record, self.coords_list)
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
        pair_list = other.coord_chop(self.length, self.size, 'exact_size')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertEquals(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)
        self.assertGreaterEqual(first_pair_len, last_pair_len)
        
    def test_chop_maxsize_bisect(self):
        pair_list = other.coord_chop(self.length, self.size, 'maxsize_bisect')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(abs(first_pair_len-last_pair_len), 1)
        
    def test_chop_maxsize_divisor(self):
        pair_list = other.coord_chop(self.length, self.size, 'maxsize_divisor')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertGreaterEqual(len(pair_list), self.length/self.size)
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)

    def test_chop_count_divisor(self):
        pair_list = other.coord_chop(self.length, self.count, 'count_divisor')
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertLessEqual(abs(len(pair_list)-self.count), 1)
        self.assertGreaterEqual(first_pair_len, last_pair_len)

    def test_chop_none(self):
        pair_list = other.coord_chop(self.length, None, None)
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        self.assertIs(len(pair_list), 1)
        self.assertEqual(first_pair_len, self.length)

    def test_recursive_bisector(self):
        pair_list = other.recursive_bisector(self.coords_list, self.size)
        print pair_list
        first_pair_len = pair_list[0][1]-pair_list[0][0]
        last_pair_len = pair_list[-1][1]-pair_list[-1][0]
        self.assertLessEqual(first_pair_len, self.size)
        self.assertLessEqual(last_pair_len, self.size)

class test_list_logic(TestCase):

    def setUp(self):
        self.list = 'spam', 'spam', 'bacon', 'spam', 'eggs', 'bacon'

    def tearDown(self):
        pass

    def test_uniqify(self):
        nuked_list = other.uniqify(self.list)
        self.assertIs(len(nuked_list), 3)

class test_os_logic(TestCase):

    def setUp(self):
        self.temp_path = 'tests/temp_dir/'

    def tearDown(self):
        try: os.rmdir(self.temp_path)
        except Exception as message: print message

    def test_ensure_dir(self):
        abs_path, report = other.ensure_dir(self.temp_path)
        existence = os.path.exists(abs_path)
        self.assertIs(report['status'], 0)
        self.assertIs(existence, True)
