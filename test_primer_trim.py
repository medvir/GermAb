#!/usr/local/bin/python3.4
'''
Script to unit test primer_trim.py. Run with nose2
'''

import os
import sys
import unittest

primer_trim_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, primer_trim_dir)
mod = __import__('primer_trim')
sys.modules["primer_trim"] = mod
import primer_trim 

class TestPrimerTrim(unittest.TestCase):
    '''test class'''
    def test_hamming_dist(self):
        '''test for hamming_dist'''
        str1, str2 = 'ACGT', 'ACGT'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)
        
        str1, str2 = 'ACGT', 'ACGg'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 1)
        
        str1, str2 = 'ACGT', 'ACtg'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 2)
        
        str1, str2 = 'ACGr', 'ACGa'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)
        
        str1, str2 = 'ACGr', 'ACGg'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)
        
        str1, str2 = 'ACGa', 'ACGr'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)
        
        str1, str2 = 'ACGg', 'ACGr'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)
        
        str1, str2 = 'AcGT', 'ACgT'
        self.assertEqual(primer_trim.hamming_dist(str1, str2), 0)


    def test_trim_fwd_primer(self):
        '''test for trim_fwd_primer'''
        seq, fwd, d = 'ACGTNNN', ['ACGT', 'TATA'], 0
        self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'NNN')
        
        seq, fwd, d = 'ACGgNNN', ['ACGT'], 1
        self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'NNN')
        
        seq, fwd, d = 'ACGTNNN', ['ACGg'], 1
        self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'NNN')
        
        seq, fwd, d = 'ACGTNNN', ['AggT'], 0
        self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'ACGTNNN')
        
        seq, fwd, d = 'ACGTNNN', ['TATA', 'ACGT'], 0
        self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'NNN')
        
        #seq, fwd, d = 'ACGTNNN', ['ACGW'], 0
        #self.assertEqual(primer_trim.trim_fwd_primer(seq, fwd, d), 'NNN')
    
        
    def test_alignment_penalty(self):
        '''test for alignment_penalty'''
        str1, str2 = 'ACGT', 'ACGT'
        self.assertEqual(primer_trim.alignment_penalty(str1, str2), 0)
        
        str1, str2 = 'ACgT', 'ACGT'
        self.assertEqual(primer_trim.alignment_penalty(str1, str2), 0)
        
        str1, str2 = 'ACGT', 'AGGT'
        self.assertEqual(primer_trim.alignment_penalty(str1, str2), 1)
        
        str1, str2 = 'ACGT', 'AGTt'
        self.assertEqual(primer_trim.alignment_penalty(str1, str2), 1)
        
        str1, str2 = 'RCGT', 'ACGT'
        self.assertEqual(primer_trim.alignment_penalty(str1, str2), 1)
        
        
    def test_count_wobbles(self):
        '''test for count_wobbles'''
        seq = 'ACGTC'
        self.assertEqual(primer_trim.count_wobbles(seq), 0)
        
        seq = 'AcGTR'
        self.assertEqual(primer_trim.count_wobbles(seq), 1)
        
        seq = 'AbCGTR'
        self.assertEqual(primer_trim.count_wobbles(seq), 2)