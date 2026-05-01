import eugene.solvers.base_min_resources_mip as emip
import unittest


class TestClass(unittest.TestCase):
    """Test case docstring."""

    def test_generate_masks(self):
        for n_loci in range(2, 7):
            results = sorted(emip.generate_masks(n_loci))
            self.assertEqual(len(results), 1 + 2**n_loci * n_loci)
            self.assertTrue(
                all(
                    results[i] != results[i + 1]
                    for i in range(len(results) - 1)
                )
            )
            mask = (1 << n_loci) - 1
            self.assertTrue(all(x_1 | x_2 == mask for x_1, x_2 in results))

            def exists_spc(x_1, x_2) -> bool:
                for k in range(n_loci):
                    mask_l = mask >> n_loci - k << n_loci - k
                    mask_r = mask >> k
                    recombination_1 = (x_1 & mask_l) | (x_2 & mask_r)
                    recombination_2 = (x_2 & mask_l) | (x_1 & mask_r)
                    if recombination_1 == mask or recombination_2 == mask:
                        return True
                print(x_1, x_2)
                return False

            self.assertTrue(all(exists_spc(x_1, x_2) for x_1, x_2 in results))

    def test_generate_recombinations(self):
        for n_loci in range(2, 6):
            results = sorted(emip.generate_recombinations(n_loci))
            self.assertEqual(
                len(results), 2**n_loci * (1 + 2**n_loci * n_loci)
            )
            self.assertTrue(
                all(
                    results[i] != results[i + 1]
                    for i in range(len(results) - 1)
                )
            )

            mask = (1 << n_loci) - 1

            def exists_spc(x_1, x_2, z_1) -> bool:
                for k in range(n_loci):
                    mask_l = mask >> n_loci - k << n_loci - k
                    mask_r = mask >> k
                    recombination_1 = (x_1 & mask_l) | (x_2 & mask_r)
                    recombination_2 = (x_2 & mask_l) | (x_1 & mask_r)
                    if recombination_1 == z_1 or recombination_2 == z_1:
                        return True
                print(x_1, x_2, z_1)
                return False

            self.assertTrue(
                all(
                    exists_spc(x >> n_loci, x & mask, z_1)
                    for x, z_1 in results
                )
            )

    def test_generate_crossings(self):
        for n_loci in range(2, 5):
            results = sorted(emip.generate_crossings(n_loci))
            mask = (1 << n_loci) - 1
            self.assertEqual(
                len(results), (2**n_loci * (1 + 2**n_loci * n_loci)) ** 2
            )
            self.assertTrue(
                all(
                    results[i] != results[i + 1]
                    for i in range(len(results) - 1)
                )
            )

            def exists_spc(x_1, x_2, z_1) -> bool:
                for k in range(n_loci):
                    mask_l = mask >> n_loci - k << n_loci - k
                    mask_r = mask >> k
                    recombination_1 = (x_1 & mask_l) | (x_2 & mask_r)
                    recombination_2 = (x_2 & mask_l) | (x_1 & mask_r)
                    if recombination_1 == z_1 or recombination_2 == z_1:
                        return True
                print(x_1, x_2, z_1)
                return False

            self.assertTrue(
                all(
                    exists_spc(x >> n_loci, x & mask, z >> n_loci)
                    and exists_spc(y >> n_loci, y & mask, z & mask)
                    for x, y, z in results
                )
            )

    def test_mip_1_locus(self):
        self.assertNotEqual(
            emip.mincross_breeding_program_distribute(dist_array=[0]),
            None,
        )

    def test_mip_2_loci(self):
        self.assertNotEqual(
            emip.mincross_breeding_program_distribute(dist_array=[0, 1]),
            None,
        )
