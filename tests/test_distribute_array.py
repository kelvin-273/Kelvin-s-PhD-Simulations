import unittest
import eugene.utils as eu
from eugene.plant_models.plant2 import PlantSPC


class TestClass(unittest.TestCase):
    def test_distribute_to_plants(self):
        self.assertEqual(
            eu.distribute_to_plants([0, 1]),
            [PlantSPC(2, 2, 2), PlantSPC(2, 1, 1)],
        )

        self.assertEqual(
            eu.distribute_to_plants([0, 1, 0]),
            [PlantSPC(3, 5, 5), PlantSPC(3, 2, 2)],
        )

        self.assertEqual(
            eu.distribute_to_plants([0, 1, 0, 2]),
            [
                PlantSPC(4, 10, 10),
                PlantSPC(4, 4, 4),
                PlantSPC(4, 1, 1),
            ],
        )
