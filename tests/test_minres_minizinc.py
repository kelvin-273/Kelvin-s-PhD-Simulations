import unittest
import eugene.solvers.base_min_resources_minizinc as ermz
import eugene.utils as eu


class TestMinRes(unittest.TestCase):

    def test_basic_is_not_none(self):
        n_loci = 2
        parents = eu.distribute_to_plants([0, 1])
        recombination_rate = [0.5]
        gamma = 0.95
        result = ermz.breeding_program(
            n_loci, parents, recombination_rate, gamma
        )
        self.assertIsNotNone(result)

    def test_basic(self):
        n_loci = 2
        parents = eu.distribute_to_plants([0, 1])
        recombination_rate = [0.5]
        gamma = 0.95
        result = ermz.breeding_program(
            n_loci, parents, recombination_rate, gamma
        )
        self.assertEqual(sum(t == "Node" for t in result.tree_type), 2)
        self.assertEqual(sum(t == "Leaf" for t in result.tree_type), 2)
