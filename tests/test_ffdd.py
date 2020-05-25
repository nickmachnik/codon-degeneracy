#!/usr/bin/env python

import unittest
import site_degeneracy


class Tests(unittest.TestCase):
    def test_site_degeneracy(self):
        res = site_degeneracy._site_degeneracy(["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"])

        self.assertTrue(all(i in res for i in [2, 0]))
        self.assertTrue(all(e in res[2] for e in [
            {'AGG', 'AGA'},
            {'CGA', 'CGU', 'CGG', 'CGC'}]))
        self.assertTrue(all(e in res[0] for e in [
            {'AGA', 'CGA'}, {'CGG', 'AGG'}]))

    def test_sub_rate(self):
        seq_a = (
            "ATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTC"
            "TAATCGCAATGGCATTCCTAATGCTTACCGAACGA")
        seq_b = (
            "ATGACCACAGTAAATCTCCTACTTATAATCATACCCACAT"
            "TAGCCGCCATAGCATTTCTCACACTCGTTGAACGA")
        subrate, _, = site_degeneracy.substitution_rate_at_ffds(
            seq_a, seq_b,
            "Vertebrate Mitochondrial", "Vertebrate Mitochondrial")
        self.assertEqual(subrate, (5, 10))


if __name__ == '__main__':
    unittest.main()
