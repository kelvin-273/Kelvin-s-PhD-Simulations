import unittest
import eugene.solvers.base_min_crossings_distribute_astar as eda
import eugene.utils as eu


class TestDistAstar(unittest.TestCase):
    """Test case docstring."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_sanitise_dist_array(self):
        def f(case):
            return eda._simplify_dist_array(case, len(case))

        self.assertEqual(f([0, 0, 1]), [0, 1])
        self.assertEqual(f([0, 1, 1]), [0, 1])

    def test_first_full_join(self):
        def f(case):
            res_raw = eda.first_full_join_raw(case)
            res = (
                eda._simplify_dist_array(res_raw[0], len(res_raw[0]))
                if res_raw is not None
                else None
            )
            return res

        self.assertEqual(f([0, 1]), [0])
        self.assertEqual(f([0, 1, 0]), None)
        self.assertEqual(f([0, 1, 2]), [0, 1])
        self.assertEqual(f([0, 1, 0, 1]), None)
        self.assertEqual(f([0, 1, 2, 1]), [0, 1, 0])
        self.assertEqual(f([0, 1, 2, 0, 3]), [0, 1, 0, 2])

    def test_branching(self):
        def f(case):
            return sorted(
                [
                    sol.xs
                    for sol in eda.branching(
                        eda.Node.from_dist_array(case),
                        eda.BRANCHING_CTX,
                    )
                ]
            )

        """
        In the current mininzinc (2024-12-13) the length of the output array is
        always the same of the length of the input array. As a consequence,
        many of the solutions that the minizinc model returns are  identical.
        For example, given the distribute array `[1, 2, 3]`, the minizinc model
        will return `[[0, 0, 2], [0, 1, 0], [0, 1, 1]]`. However, `[0, 0, 2]`
        is structurally identical to `[0, 1, 1]` in that the any optimal
        solution to `[0, 0, 2]` exhibits the same structure to some opitmal
        solution to `[0, 1, 1]`.

        This necessitates a second layer of symmetry breaking,
        left as an exercise to a future computer scientist.

        Wait (2024-12-13), this symmetry breaking is meant to be performed when
        the `simplify_results` option is set to `True`.
        See `test_branching_simplified` below.
        """
        # TODO: Break symmetry on equivalent output arrays <13-12-24> #
        self.assertEqual(f([0, 1]), [[0]])
        self.assertEqual(f([0, 1, 0]), [[0, 1]])
        self.assertEqual(f([0, 1, 2]), [[0, 1]])
        self.assertEqual(f([0, 1, 0, 1]), [[0, 1], [0, 1, 0]])
        self.assertEqual(f([0, 1, 0, 2]), [[0, 1, 0]])
        self.assertEqual(f([0, 1, 2, 0]), [[0, 1, 0]])
        self.assertEqual(f([0, 1, 2, 1]), [[0, 1, 0]])
        self.assertEqual(f([0, 1, 2, 3]), [[0, 1, 2]])
        self.assertEqual(
            f([0, 1, 0, 1, 0]), [[0, 1, 0, 1], [0, 1, 2], [0, 1, 2, 0]]
        )
        self.assertEqual(f([0, 1, 0, 1, 2]), [[0, 1, 0, 1]])
        self.assertEqual(
            # Verified that these solutions are consistent but not total
            f([0, 1, 0, 2, 0]),
            [
                [0, 1, 0, 1, 0],
                [0, 1, 0, 2],
                [0, 1, 0, 2, 1],
                [0, 1, 2, 0],
                [0, 1, 2, 0, 2],
                [0, 1, 2, 1],
            ],
        )
        self.assertEqual(f([0, 1, 0, 2, 1]), [[0, 1, 0, 1]])
        self.assertEqual(f([0, 1, 0, 2, 3]), [[0, 1, 0, 2]])
        self.assertEqual(
            f([0, 1, 2, 0, 1]),
            [
                [0, 1, 0, 2],
                [0, 1, 0, 2, 1],
                [0, 1, 2],
                [0, 1, 2, 0],
                [0, 1, 2, 0, 2],
                [0, 1, 2, 1],
                [0, 1, 2, 1, 0],
            ],
        )
        self.assertEqual(f([0, 1, 2, 0, 2]), [[0, 1, 0, 1]])
        self.assertEqual(f([0, 1, 2, 0, 3]), [[0, 1, 0, 2]])
        self.assertEqual(
            f([0, 1, 2, 1, 0]),
            [
                [0, 1, 0, 1, 2],
                [0, 1, 0, 2],
                [0, 1, 2],
                [0, 1, 2, 0],
                [0, 1, 2, 0, 1],
                [0, 1, 2, 1],
                [0, 1, 2, 1, 2],
            ],
        )
        self.assertEqual(f([0, 1, 2, 1, 2]), [[0, 1, 0, 1]])
        self.assertEqual(f([0, 1, 2, 1, 3]), [[0, 1, 0, 2]])
        self.assertEqual(f([0, 1, 2, 3, 0]), [[0, 1, 2, 0]])
        self.assertEqual(f([0, 1, 2, 3, 1]), [[0, 1, 2, 0]])
        self.assertEqual(f([0, 1, 2, 3, 2]), [[0, 1, 2, 1]])
        self.assertEqual(f([0, 1, 2, 3, 4]), [[0, 1, 2, 3]])

        self.assertEqual(
            f([0, 1, 2, 0, 1]),
            sorted(
                [
                    [0, 1, 2, 0],
                    [0, 1, 2],
                    [0, 1, 0, 2, 1],
                    [0, 1, 0, 2],
                    [0, 1, 2, 0, 2],
                    [0, 1, 2, 1, 0],
                    [0, 1, 2, 1],
                ]
            ),
        )

        self.assertTrue([0, 1, 2, 2, 0, 1] not in f([0, 1, 0, 2, 0, 3, 0, 1]))

    def test_generate_raw_redistributions(self):
        instance = [0, 1, 2, 0, 1]
        raw_redistributions = eda.generate_raw_redistributions(instance)
        # self.assertTrue(
        #     eda.is_redistribution_of([0, 1, 2, 1, 0], [0, 1, 2, 0, 1])
        # )
        self.assertTrue([0, 1, 2, 1, 0] in raw_redistributions)

    def test_branching_simplified_always_distribute(self):
        for n_loci in range(6, 11):
            for _ in range(100):
                xs = eu.random_distribute_instance(n_loci)
                children = eda.branching(
                    eda.Node.from_dist_array(xs),
                    eda.BRANCHING_CTX,
                )
            for child in children:
                self.assertTrue(
                    eda.is_distribute(child.xs),
                    msg=f"child {child.xs} of parent {xs} is non-distribute",
                )

    def test_astar_objective(self):
        cases = [
            {"check": 4, "xs": [0, 1, 2, 0, 1]},
            {"check": 4, "xs": [0, 1, 2, 1, 0]},
            {"check": 4, "xs": [0, 1, 0, 1, 0, 1]},
            {"check": 5, "xs": [0, 1, 2, 0, 1, 3]},
            {"check": 5, "xs": [0, 1, 2, 0, 2, 1]},
            {"check": 5, "xs": [0, 1, 2, 1, 0, 3]},
            {"check": 5, "xs": [0, 1, 2, 1, 3, 0]},
            {"check": 5, "xs": [0, 1, 2, 3, 0, 1]},
            {"check": 5, "xs": [0, 1, 2, 3, 1, 0]},
            {"check": 5, "xs": [0, 1, 2, 3, 1, 2]},
            {"check": 5, "xs": [0, 1, 2, 3, 2, 0]},
            {"check": 5, "xs": [0, 1, 2, 3, 2, 1]},
            {"check": 5, "xs": [0, 1, 0, 1, 0, 2, 1]},
            {"check": 5, "xs": [0, 1, 0, 1, 2, 0, 2]},
            {"check": 5, "xs": [0, 1, 0, 2, 0, 1, 2]},
            {"check": 5, "xs": [0, 1, 0, 2, 1, 2, 0]},
            {"check": 5, "xs": [0, 1, 2, 0, 2, 0, 2]},
            {"check": 5, "xs": [0, 1, 2, 0, 2, 1, 2]},
            {"check": 5, "xs": [0, 1, 2, 1, 0, 2, 0]},
            {"check": 5, "xs": [0, 1, 2, 1, 2, 1, 0]},
            {"check": 5, "xs": [0, 1, 2, 1, 2, 1, 2]},
            {"check": 6, "xs": [0, 1, 0, 2, 3, 0, 2]},
            {"check": 6, "xs": [0, 1, 0, 2, 3, 2, 0]},
            {"check": 6, "xs": [0, 1, 2, 0, 1, 3, 4]},
            {"check": 6, "xs": [0, 1, 2, 0, 2, 3, 1]},
            {"check": 6, "xs": [0, 1, 2, 0, 3, 1, 0]},
            {"check": 6, "xs": [0, 1, 2, 1, 0, 3, 0]},
            {"check": 6, "xs": [0, 1, 2, 1, 3, 1, 0]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 2, 1]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 3, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 4, 1]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 4, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 1, 2, 4]},
            {"check": 6, "xs": [0, 1, 2, 3, 1, 4, 0]},
            {"check": 6, "xs": [0, 1, 2, 3, 1, 4, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 2, 0, 4]},
            {"check": 6, "xs": [0, 1, 2, 3, 2, 4, 1]},
            {"check": 6, "xs": [0, 1, 2, 3, 4, 0, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 4, 0, 3]},
            {"check": 6, "xs": [0, 1, 2, 3, 4, 2, 3]},
            {"check": 6, "xs": [0, 1, 2, 3, 4, 3, 0]},
            {"check": 6, "xs": [0, 1, 2, 3, 4, 3, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 1, 4, 1]},
            {"check": 6, "xs": [0, 1, 0, 2, 3, 0, 1, 0]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 4, 0, 5]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 5, 2, 1]},
            {"check": 6, "xs": [0, 1, 0, 2, 3, 1, 2, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 3, 5, 2]},
            {"check": 6, "xs": [0, 1, 2, 1, 3, 2, 3, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 3, 0, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 0, 4, 1]},
            {"check": 6, "xs": [0, 1, 0, 2, 3, 2, 0, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 1, 0, 2, 1]},
            {"check": 6, "xs": [0, 1, 2, 3, 2, 0, 1, 3]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 4, 2, 1]},
            {"check": 6, "xs": [0, 1, 2, 0, 2, 3, 1, 2]},
            {"check": 6, "xs": [0, 1, 2, 1, 3, 2, 3, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 5, 1, 2]},
            {"check": 6, "xs": [0, 1, 2, 0, 3, 1, 0, 3]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 0, 1, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 2, 4, 5]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 0, 1, 5]},
            {"check": 6, "xs": [0, 1, 2, 1, 2, 0, 3, 1]},
            {"check": 6, "xs": [0, 1, 0, 1, 2, 1, 3, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 1, 4, 5]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 0, 4, 1]},
            {"check": 6, "xs": [0, 1, 2, 0, 2, 3, 2, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 3, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 2, 4, 3]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 5, 3, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 5, 2, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 4, 0, 1]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 2, 0, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 0, 3, 0, 4]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 2, 0, 1]},
            {"check": 6, "xs": [0, 1, 2, 0, 1, 2, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 1, 4, 0]},
            {"check": 6, "xs": [0, 1, 0, 1, 2, 3, 0, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 1, 2, 3, 2]},
            {"check": 6, "xs": [0, 1, 2, 3, 2, 3, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 0, 3, 4, 0, 3]},
            {"check": 6, "xs": [0, 1, 2, 3, 0, 1, 0, 3]},
            {"check": 7, "xs": [0, 1, 0, 2, 0, 3, 4, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 3, 2, 0]},
            {"check": 6, "xs": [0, 1, 2, 1, 0, 2, 1, 2]},
            {"check": 6, "xs": [0, 1, 2, 1, 2, 3, 1, 0]},
            {"check": 6, "xs": [0, 1, 2, 1, 2, 3, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 3, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 1, 4, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 4, 0, 4, 3]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 3, 4, 3, 1]},
            {"check": 7, "xs": [0, 1, 0, 1, 2, 0, 2, 1, 0]},
            {"check": 7, "xs": [0, 1, 2, 0, 2, 3, 1, 2, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 1, 2, 1, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 4, 1, 3, 1]},
            {"check": 7, "xs": [0, 1, 0, 2, 0, 1, 3, 1, 3]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 3, 2, 0, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 1, 2, 0, 3]},
            {"check": 8, "xs": [0, 1, 2, 3, 0, 4, 5, 1, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 2, 0, 1, 3]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 0, 1, 0, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 2, 0, 1, 2]},
            {"check": 8, "xs": [0, 1, 2, 0, 3, 1, 4, 0, 5]},
            {"check": 7, "xs": [0, 1, 0, 1, 2, 3, 0, 2, 3]},
            {"check": 7, "xs": [0, 1, 0, 2, 1, 3, 4, 3, 2]},
            {"check": 8, "xs": [0, 1, 2, 0, 3, 4, 5, 3, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 1, 2, 4, 0]},
            {"check": 8, "xs": [0, 1, 2, 3, 0, 4, 5, 4, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 1, 3, 4, 2]},
            {"check": 7, "xs": [0, 1, 2, 1, 0, 3, 0, 2, 3]},
            {"check": 7, "xs": [0, 1, 2, 3, 1, 0, 4, 2, 3]},
            {"check": 7, "xs": [0, 1, 0, 2, 3, 2, 1, 3, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 2, 0, 1, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 1, 0, 2, 4]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 3, 1, 2, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 1, 3, 0, 1, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 1, 4, 0, 1]},
            {"check": 7, "xs": [0, 1, 0, 2, 3, 1, 4, 1, 3]},
            {"check": 8, "xs": [0, 1, 2, 3, 1, 4, 5, 0, 6]},
            {"check": 7, "xs": [0, 1, 2, 0, 3, 2, 1, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 2, 3, 1, 2]},
            {"check": 8, "xs": [0, 1, 2, 3, 4, 5, 6, 5, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 2, 1, 2, 0]},
            {"check": 6, "xs": [0, 1, 0, 2, 0, 1, 0, 2, 1]},
            {"check": 7, "xs": [0, 1, 2, 1, 2, 3, 1, 4, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 1, 3, 1, 2, 0]},
            {"check": 7, "xs": [0, 1, 0, 1, 2, 1, 3, 0, 4]},
            {"check": 7, "xs": [0, 1, 2, 1, 0, 3, 0, 4, 3]},
            {"check": 7, "xs": [0, 1, 2, 0, 3, 2, 1, 2, 0]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 2, 1, 0, 4]},
            {"check": 7, "xs": [0, 1, 0, 1, 2, 1, 3, 1, 4]},
            {"check": 7, "xs": [0, 1, 0, 2, 0, 2, 3, 1, 2]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 4, 3, 0, 4]},
            {"check": 8, "xs": [0, 1, 2, 3, 4, 3, 5, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 2, 4, 1, 4]},
            {"check": 7, "xs": [0, 1, 2, 0, 3, 2, 4, 3, 2]},
            {"check": 8, "xs": [0, 1, 2, 0, 3, 0, 4, 0, 5]},
            {"check": 7, "xs": [0, 1, 0, 2, 0, 3, 1, 4, 3]},
            {"check": 6, "xs": [0, 1, 0, 2, 1, 2, 1, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 0, 3, 1, 4]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 0, 3, 2, 1]},
            {"check": 8, "xs": [0, 1, 2, 3, 1, 4, 1, 5, 1]},
            {"check": 8, "xs": [0, 1, 0, 2, 3, 2, 4, 5, 4]},
            {"check": 8, "xs": [0, 1, 2, 3, 1, 4, 5, 0, 1]},
            {"check": 7, "xs": [0, 1, 2, 0, 2, 3, 4, 2, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 1, 3, 2, 0]},
            {"check": 7, "xs": [0, 1, 2, 3, 0, 4, 3, 0, 2]},
            {"check": 7, "xs": [0, 1, 2, 0, 3, 1, 4, 3, 2]},
            {"check": 7, "xs": [0, 1, 2, 1, 3, 1, 3, 4, 0]},
            {"check": 7, "xs": [0, 1, 2, 1, 2, 3, 0, 2, 4]},
            {"check": 7, "xs": [0, 1, 2, 3, 4, 1, 0, 4, 2]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 3, 0, 2, 4]},
            {"check": 7, "xs": [0, 1, 0, 2, 3, 4, 1, 3, 1]},
            {"check": 7, "xs": [0, 1, 2, 3, 2, 4, 2, 1, 3]},
            {"check": 7, "xs": [0, 1, 2, 0, 1, 2, 3, 4, 1]},
        ]
        for case in cases:
            self.assertEqual(
                eda.breeding_program_distribute(case["xs"]).objective,
                case["check"],
                msg=f"failed for {case['xs']}",
            )


if __name__ == "__main__":
    zss = eda.generate_simplified_redistributions_brute_force(
        [0, 1, 2, 1, 0, 3, 4, 2, 3, 2]
    )
    for zs in zss:
        print(
            f"segs: {len(zs)},\tgams: {max(zs) + 1},\ttotal: {len(zs) + max(zs) + 1},\tzs: {zs}"
        )
    __import__("pprint").pprint(sorted(zss))
