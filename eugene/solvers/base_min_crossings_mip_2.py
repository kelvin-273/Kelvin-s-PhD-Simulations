from typing import List
from eugene.solution import BaseSolution
from eugene.plant_models.plant2 import PlantSPC
import gurobipy as gp


def breeding_program_distribute(dist_array: List[int]) -> BaseSolution:
    n_loci = len(dist_array)
    pop_0 = None
    return breeding_program(n_loci, pop_0)


def breeding_program(n_loci: int, pop_0: List[PlantSPC]) -> BaseSolution:
    n_pop = len(pop_0)
    m = gp.Model()
    T = upper_bound(n_loci, pop_0)
    tree_data = [
        [m.addVar(ub=(1 << n_loci) - 1) for _ in range(2)]
        for _ in range(T + n_pop)
    ]

    # first |pop_0| cells of tree_data are allocated to pop_0

    # each gamete of each crossing must be created by recombination
    # of some previous genotype

    # count the number of actual crossings

    m.optimize()

    return BaseSolution(
        tree_data=[],
        tree_type=[],
        tree_left=[],
        tree_right=[],
        objective=m.objVal,
    )


def upper_bound(n_loci: int, pop_0: List[PlantSPC]) -> int:
    return n_loci


if __name__ == "__main__":
    n_loci = 2
    print(
        breeding_program(
            n_loci,
            [PlantSPC(n_loci, 0b10, 0b10,), PlantSPC(n_loci, 0b01, 0b01,),],
        )
    )
