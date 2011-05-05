__author__ = 'GG'

import sys
sys.path.append("/Users/GG/codespace/trappist")

import os
from unittest import TestCase
from Bio import AlignIO
from Bio.Seq import Seq, SeqRecord
from analysis import sequence_file_ops, alignment

class test_alignments(TestCase):

    def setUp(self):
        # create temp directory
        self.temp_dir = 'tests/temp_data/'
        os.mkdir(self.temp_dir)
        # create some sequence records
        self.seq_1 = Seq('AATTTAATGCATATCTAGGAGAGAGATTTTTGGCGCTAACGGGG')
        self.seq_2 = Seq('GGTTAATGGCGCAGGCTAGGAGAGAGATTTTTGGCGCTCGAAGCGGACG')
        self.seq_3 = Seq('GGTTAATGCCGCAGGCTAGGAGAGAGATGTGGCGCTCGCGGCGGCG')
        self.record_1 = SeqRecord(self.seq_1, id='temp_1')
        self.record_2 = SeqRecord(self.seq_2, id='temp_2')
        self.record_3 = SeqRecord(self.seq_3, id='temp_3')
        # create record sets
        self.two_records = [self.record_1, self.record_2]
        self.three_records = [self.record_1, self.record_2, self.record_3]
        # define file names
        self.seq_file_2 = self.temp_dir+'temp_2_in.fas'
        self.seq_file_3 = self.temp_dir+'temp_3_in.fas'
        self.aln_file_2 = self.temp_dir+'temp_2_in.aln'
        self.aln_file_3 = self.temp_dir+'temp_3_in.aln'
        self.dnd_file_2 = self.temp_dir+'temp_2_in.dnd'
        self.dnd_file_3 = self.temp_dir+'temp_3_in.dnd'
        self.muscle_log = self.temp_dir+'temp_log.txt'
        self.seq_file_m1 = self.temp_dir+'temp_m1_in.fas'
        self.seq_file_m2 = self.temp_dir+'temp_m2_in.fas'
        self.seq_file_m3 = self.temp_dir+'temp_m3_in.fas'
        self.mauve_out = self.temp_dir+'temp_out.mauve'
        self.mauve_out_bb = self.mauve_out+'.backbone'
        self.mauve_out_bbcols = self.mauve_out+'.bbcols'
        self.mauve_sslist_1 = self.seq_file_m1+'.sslist'
        self.mauve_sslist_2 = self.seq_file_m2+'.sslist'
        self.mauve_sslist_3 = self.seq_file_m3+'.sslist'
        # create temp files
        sequence_file_ops.write_fasta(self.seq_file_2, self.two_records)
        sequence_file_ops.write_fasta(self.seq_file_3, self.three_records)
        sequence_file_ops.write_fasta(self.seq_file_m1, self.record_1)
        sequence_file_ops.write_fasta(self.seq_file_m2, self.record_2)
        sequence_file_ops.write_fasta(self.seq_file_m3, self.record_3)
        # command line for Mauve
        self.cline = "progressiveMauve --output="+self.mauve_out \
                     +" "+self.seq_file_m1+" "+self.seq_file_m2 \
                     +" "+self.seq_file_m3

    def tearDown(self):
        temp_files = [self.seq_file_2, self.seq_file_3, self.aln_file_2,
                      self.aln_file_3, self.dnd_file_2, self.dnd_file_3,
                      self.muscle_log, self.seq_file_m1, self.seq_file_m2,
                      self.seq_file_m3, self.mauve_out, self.mauve_out_bb,
                      self.mauve_out_bbcols, self.mauve_sslist_1,
                      self.mauve_sslist_2, self.mauve_sslist_3]
        for file in temp_files:
            try: os.remove(file)
            except Exception as message: print message
        try: os.rmdir(self.temp_dir)
        except Exception as message: print message

    def test_align_clustal_2(self):
        # prepare alignment
        clustal_report = alignment.align_clustal(self.seq_file_2)
        self.assertIsNot(clustal_report['output'], '')
        self.assertIsNone(clustal_report['error'])
        # test aligned output
        aligned = AlignIO.read(open(self.aln_file_2), 'clustal')
        aln_length = aligned.get_alignment_length()
        self.assertIs(aln_length, 49)

    def test_align_clustal_3(self):
        # prepare alignment
        clustal_report = alignment.align_clustal(self.seq_file_3)
        self.assertIsNot(clustal_report['output'], '')
        self.assertIsNone(clustal_report['error'])
        # test aligned output
        aligned = AlignIO.read(open(self.aln_file_3), 'clustal')
        aln_length = aligned.get_alignment_length()
        self.assertIs(aln_length, 50)

    def test_parse_clustal_idstars_2(self):
        alignment.align_clustal(self.seq_file_2)
        idntot = alignment.parse_clustal_idstars(self.aln_file_2)
        self.assertIs(idntot, 31)

    def test_parse_clustal_idstars_3(self):
        alignment.align_clustal(self.seq_file_3)
        idntot = alignment.parse_clustal_idstars(self.aln_file_3)
        self.assertIs(idntot, 30)

    def test_align_muscle_2(self):
        # prepare alignment
        muscle_report = alignment.align_muscle(self.seq_file_2,
                                               self.aln_file_2,
                                               self.muscle_log)
        self.assertIsNot(muscle_report['output'], '')
        self.assertIsNone(muscle_report['error'])
        # test aligned output
        aligned = AlignIO.read(open(self.aln_file_2), 'clustal')
        aln_length = aligned.get_alignment_length()
        self.assertIs(aln_length, 49)

    def test_align_muscle_3(self):
        # prepare alignment
        muscle_report = alignment.align_muscle(self.seq_file_3,
                                               self.aln_file_3,
                                               self.muscle_log)
        self.assertIsNot(muscle_report['output'], '')
        self.assertIsNone(muscle_report['error'])
        # test aligned output
        aligned = AlignIO.read(open(self.aln_file_3), 'clustal')
        aln_length = aligned.get_alignment_length()
        self.assertIs(aln_length, 49)

    def test_align_mauve(self):
        mauve_report = alignment.align_mauve(self.cline)
        self.assertIsNot(mauve_report['output'], '')
        self.assertIsNone(mauve_report['error'])
        existence = os.path.exists(self.mauve_out_bb)
        self.assertIs(existence, True)
