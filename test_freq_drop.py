#!/usr/local/bin/python3.4
'''
Script to unit test freq_drop.py. Run with nose2
'''

import os
import sys
import unittest

freq_drop_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, freq_drop_dir)
mod = __import__('freq_drop')
sys.modules["freq_drop"] = mod
import freq_drop

class TestFreqDrop(unittest.TestCase):
    '''test class'''
    def test_freq_drop(self):
        '''test for freq_drop'''
        reads = '1000', '99'
        self.assertEqual(freq_drop.freq_drop(reads), 1)
        
        reads = '1000', '999'
        self.assertEqual(freq_drop.freq_drop(reads), 2)
        
        reads = '1000', '999', '990'
        self.assertEqual(freq_drop.freq_drop(reads), 3)

        reads = '1000', '990', '9'
        self.assertEqual(freq_drop.freq_drop(reads), 2)

        reads = '1000', '990', '90', '88', '77'
        self.assertEqual(freq_drop.freq_drop(reads), 2)
        
        reads = '1076', '1006', '291', '207', '121', '28' , '25'
        self.assertEqual(freq_drop.freq_drop(reads), 2)