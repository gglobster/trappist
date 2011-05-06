__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from analysis import text_manipulation

class test_text_manipulation(TestCase):

    def setUp(self):
        # create temp directory
        self.temp_dir = "tests/temp_data/"
        os.mkdir(self.temp_dir)
        # create content for the temporary file
        self.header_contents = ['A.these', 'B.are', 'C.column', 'D.headers']
        self.line_1_contents = ['item_1_A', 'item_1_B', 'item_1_C', 'item_1_D']
        self.line_2_contents = ['item_2_A', 'item_2_B', 'item_2_C', 'item_2_D']
        self.line_3_contents = ['item_3_A', 'item_3_B', 'item_3_C', 'item_3_D']
        self.raw_contents = [self.header_contents, self.line_1_contents,
                             self.line_2_contents, self.line_3_contents]
        # transform into a single string
        self.str_contents = ""
        for line_set in self.raw_contents:
            joined = "\t".join(line_set)
            self.str_contents += joined+"\n"
        # write to file
        self.filename = self.temp_dir+"temp_file.txt"
        self.temp_file = open(self.filename, 'w')
        self.temp_file.write(self.str_contents)
        self.temp_file.close()

    def tearDown(self):
        try: os.remove(self.filename)
        except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_td_txt_file_load_no_skip(self):
        rawlines = text_manipulation.td_txt_file_load(self.filename, 0)
        self.assertEqual(len(rawlines), 4)
        rejoined = ''.join(rawlines)
        self.assertEqual(rejoined, self.str_contents)

    def test_td_txt_file_load_skip_1(self):
        rawlines = text_manipulation.td_txt_file_load(self.filename, 1)
        self.assertEqual(len(rawlines), 3)
        rejoined = ''.join(rawlines)
        self.assertNotEqual(rejoined, self.str_contents)

    def test_strippa_clean(self):
        rawlines = text_manipulation.td_txt_file_load(self.filename, 0)
        line_tuple = text_manipulation.strippa_clean(rawlines[1])
        for index in range(0, 4):
            self.assertEqual(line_tuple[index], self.line_1_contents[index])

    def test_adaptive_list_load(self):
        trimlines = text_manipulation.adaptive_list_load(self.filename, 0,
                                                         (1,3))
        self.assertEqual(len(trimlines), 4)
        self.assertEqual(trimlines[1][1], self.line_1_contents[3])
