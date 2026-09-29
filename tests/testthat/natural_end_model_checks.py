"""Behaviour checks for the local heterogeneous fragment generator."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/benchmarks/rust'))
from natural_end_model import *

class NaturalEndModelChecks(unittest.TestCase):
    def test_geometry_and_reference(self):
        seq='ACGTTGCA'*40
        sites=np.array([60,90,120])
        c=fragment_candidates(seq,sites,np.array([1.,2.,3.]))
        np.testing.assert_array_equal(c['start']+c['offset'],c['site'])
        self.assertEqual(set(c['length']),set(LENGTHS))
        for start,width,read in zip(c['start'],c['length'],c['sequence']):
            self.assertEqual(seq[start:start+width],read)
        for length,prob,offset in zip(LENGTHS,LENGTH_PROB,MODAL_OFFSETS):
            selected=c['length']==length
            self.assertAlmostEqual(c['weight'][selected].sum()/c['weight'].sum(),prob)
            self.assertEqual(set(c['offset'][selected]),set(offset+JITTER))
        with self.assertRaises(ValueError):fragment_candidates('ACGT',np.array([2]),np.array([1.]))

    def test_noise_is_reproducible_and_nonuniform(self):
        w=noisy_site_weights(np.ones(1000),np.random.default_rng(7))
        np.testing.assert_array_equal(w,noisy_site_weights(np.ones(1000),np.random.default_rng(7)))
        self.assertTrue(np.all(np.isfinite(w)&(w>0)))
        self.assertGreater(w.std(),.2)
        self.assertGreater(np.sum(w<.15),20)
        with self.assertRaises(ValueError):noisy_site_weights([0],np.random.default_rng(1))

    def test_library_factors_multiply_and_preserve_depth(self):
        c=dict(weight=np.ones(4),five_weight=np.array([2.,2.,1.,1.]),
               three_weight=np.array([2.,1.,2.,1.]))
        original=c['weight'].copy()
        cases={'control':[1,1,1,1],'five_prime':[2,2,1,1],
               'three_prime':[2,1,2,1],'both_ends':[4,2,2,1]}
        for scenario,relative in cases.items():
            sampled=sample_library(c,100000,scenario,np.random.default_rng(3))
            self.assertEqual(sampled.sum(),100000)
            np.testing.assert_allclose(sampled/100000,np.array(relative)/sum(relative),atol=.005)
        np.testing.assert_array_equal(c['weight'],original)
        with self.assertRaises(ValueError):sample_library(c,-1,'control',np.random.default_rng())

    def test_biological_dinucleotide_ends(self):
        self.assertAlmostEqual(dinucleotide_weights(['GC'],'five')[0],1.35*1.10*1.12)
        self.assertAlmostEqual(dinucleotide_weights(['TC'],'three')[0],.90*1.40*1.20)
        with self.assertRaises(ValueError):dinucleotide_weights(['NC'],'three')

if __name__=='__main__':unittest.main()
