__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

from unittest import TestCase
from analysis import dataset_sort

class test_reorder_matches(TestCase):

    def setUp(self):
        # set up basic matches of query (q) to contig (c)
        q1_c1 = {'contig_id': 'test_contig_1', 'details': {'match_p100': 80}}
        q2_c1 = {'contig_id': 'test_contig_1', 'details': {'match_p100': 30}}
        q2_c2 = {'contig_id': 'test_contig_2', 'details': {'match_p100': 60}}
        q3_c2 = {'contig_id': 'test_contig_2', 'details': {'match_p100': 90}}
        q4_c3 = {'contig_id': 'test_contig_3', 'details': {'match_p100': 70}}
        # set up match sets
        self.q1_matches = {'query_id': 'test_query_1', 'matches': [q1_c1]}
        self.q2_matches = {'query_id': 'test_query_2', 'matches': [q2_c1,
                                                                   q2_c2]}
        self.q3_matches = {'query_id': 'test_query_3', 'matches': [q3_c2]}
        self.q4_matches = {'query_id': 'test_query_4', 'matches': [q4_c3]}

    def tearDown(self):
        pass

    def test_reorder_matches_by_contig(self):
        q_matches = [self.q1_matches, self.q2_matches, self.q3_matches,
                     self.q4_matches]
        cc_matches = dataset_sort.reorder_matches_by_contig(q_matches)
        self.assertIs(len(cc_matches), 3)
        self.assertIs(len(cc_matches['test_contig_1']), 2)
        self.assertIs(len(cc_matches['test_contig_2']), 2)
        self.assertIs(len(cc_matches['test_contig_3']), 1)
        self.assertIs(cc_matches['test_contig_3'][0]['details']['match_p100'],
                      70)

