__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from analysis import other_logic

class test_list_logic(TestCase):

    def setUp(self):
        self.list = 'spam', 'spam', 'bacon', 'spam', 'eggs', 'bacon'

    def tearDown(self):
        pass

    def test_uniqify(self):
        nuked_list = other_logic.uniqify(self.list)
        self.assertIs(len(nuked_list), 3)

class test_os_logic(TestCase):

    def setUp(self):
        self.temp_path = "tests/temp_dir/"

    def tearDown(self):
        try: os.rmdir(self.temp_path)
        except Exception as message: print message

    def test_ensure_dir(self):
        abs_path, report = other_logic.ensure_dir(self.temp_path)
        existence = os.path.exists(abs_path)
        self.assertIs(report['status'], 0)
        self.assertIs(existence, True)

# TODO: make test for create_id