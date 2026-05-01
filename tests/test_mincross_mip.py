import unittest
import eugene.solvers.base_min_resources_mip as emip
import itertools


class TestMIP(unittest.TestCase):

    def test_generate_crossings_for_offspirng(self):
        n_loci = 2
        z = 15
        self.assertEqual(
            sorted(emip.generate_crossings_for_offpring(n_loci, z)),
            list(itertools.product([6, 7, 9, 11, 13, 15], repeat=2)),
        )
