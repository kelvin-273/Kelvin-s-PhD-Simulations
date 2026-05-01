import unittest
import eugene.solvers.base_min_crossings_minizinc as emzn
import eugene.plant_models.plant2 as ep2


class TestMinCrossingsMinizinc(unittest.TestCase):
    def test_mincross_minizinc_trivial(self):
        n_loci = 1
        pop_0 = ep2.PlantSPC.initial_pop_singles_homo(n_loci)
        result = emzn.breeding_program(n_loci, pop_0)
        self.assertIsNotNone(result, "result is None")

    def test_mincross_minizinc_2_loci(self):
        n_loci = 2
        pop_0 = ep2.PlantSPC.initial_pop_singles_homo(n_loci)
        result = emzn.breeding_program(n_loci, pop_0)
        self.assertIsNotNone(result, "result is None")
        self.assertEqual(result.objective, 2, "result should have 2 crossings")
