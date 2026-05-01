import contextlib
import minizinc
from math import ceil
from typing import List
from datetime import timedelta

import eugene.utils as eu
import eugene.plant_models.plant2 as ep2
from eugene.solution import BaseSolution


def breeding_program(n_loci, pop_0, recombination_rate, gamma=0.95, ctx=None, timeout=None) -> BaseSolution:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """
    items = instance_array_genotype(pop_0, max_crossovers=1).items()

    if ctx is None:
        ctx = MinizincContext.from_solver_and_model_file(
            "sat", "./eugene/solvers/minizinc/minres.mzn",
        )
    if timeout is not None:
        timeout = timedelta(seconds=timeout)
    with ctx.instance.branch() as instance:

        for k, v in items:
            instance[k] = v
        instance["recombinationRate"] = recombination_rate
        instance["gamma"] = gamma

        result = instance.solve(free_search=True, timeout=timeout)
        return BaseSolution(
            tree_data=result["treeData"],
            tree_type=result["treeType"],
            tree_left=result["treeLeft"],
            tree_right=result["treeRight"],
            objective=result.objective,
        )


class MinizincContext:
    """
    A container for the model, solver, and instance configurations for
    minizinc models.
    """

    def __init__(
        self,
        model: minizinc.Model,
        solver: minizinc.Solver,
        instance: minizinc.Instance,
    ):
        self.model = model
        self.solver = solver
        self.instance = instance

    @staticmethod
    def from_solver_and_model_file(solver: str, model_file: str):
        model = minizinc.Model(model_file)
        solver = minizinc.Solver.lookup(solver)
        instance = minizinc.Instance(solver, model)
        return MinizincContext(model, solver, instance)

    @contextlib.contextmanager
    def branch(self):
        with self.instance.branch() as child:
            yield MinizincContext(self.model, self.solver, child)


def instance_array_genotype_homo(instance_array, max_crossovers):
    """
    Creates an instance as a dictionary of an instance of homozygous genotypes
    from an instance array and an integer representing the maximum allowed
    crossovers.
    """
    n_loci = len(instance_array)
    n_gametes = max(instance_array) + 1
    return {
        "maxCrossovers": max_crossovers,
        "nLoci": n_loci,
        "nParents": n_gametes,
        "nTreeCells": n_loci + n_gametes,
        "parents": [
            [[1 if j == i else 0 for j in instance_array]] * 2
            for i in range(n_gametes)
        ],
    }


def instance_array_genotype(pop_0: List[ep2.PlantSPC], max_crossovers):
    """
    Creates an instance as a dictionary of an instance of homozygous genotypes
    from an instance array and an integer representing the maximum allowed
    crossovers.
    """
    if len(pop_0) == 0:
        raise NotImplementedError("Need a way to account for empty instances")

    n_loci = pop_0[0].n_loci
    n_pop = len(pop_0)
    return {
        "maxCrossovers": max_crossovers,
        "nLoci": n_loci,
        "nParents": n_pop,
        "nTreeCells": n_loci + n_pop,
        "parents": [x.to_bitlist() for x in pop_0],
    }


DEFAULT_CTX = MinizincContext.from_solver_and_model_file(
    "sat", "./eugene/solvers/minizinc/minres.mzn"
)


def check_solution(solution: BaseSolution) -> bool:
    n_loci = solution.n_loci
    n_plants = solution.n_plants
    n_tcells = len(solution.tree_data)
    n_cross = solution.crossings()

    def crossable(i: int) -> bool:
        return crossable_gam(
            solution.tree_data[solution.tree_left[i] - 1][0],
            solution.tree_data[solution.tree_left[i] - 1][1],
            solution.tree_data[i][0],
        ) and crossable_gam(
            solution.tree_data[solution.tree_right[i] - 1][0],
            solution.tree_data[solution.tree_right[i] - 1][1],
            solution.tree_data[i][1],
        )

    def crossable_gam(gx, gy, gz) -> bool:
        mat_pref = 0
        mat_suff = 0
        while mat_pref < n_loci and gz[mat_pref] == gx[mat_pref]:
            mat_pref += 1
        while (
            mat_suff < n_loci - mat_pref
            and gz[-1 - mat_suff] == gy[-1 - mat_suff]
        ):
            mat_suff += 1

        if mat_pref + mat_suff >= n_loci:
            return True

        mat_pref = 0
        mat_suff = 0
        while mat_pref < n_loci and gz[mat_pref] == gy[mat_pref]:
            mat_pref += 1
        while (
            mat_suff < n_loci - mat_pref
            and gz[-1 - mat_suff] == gx[-1 - mat_suff]
        ):
            mat_suff += 1

        return mat_pref + mat_suff >= n_loci

    def is_ideotype(i: int) -> bool:
        return solution.tree_data[i] == [[1] * n_loci] * 2

    assert len(solution.tree_data) == n_tcells
    assert len(solution.tree_type) == n_tcells
    assert len(solution.tree_left) == n_tcells
    assert len(solution.tree_right) == n_tcells
    # Nulls at the back
    assert n_plants <= n_tcells
    assert all(ty in ("Node", "Leaf") for ty in solution.tree_type[:n_plants])
    assert all(ty == "Null" for ty in solution.tree_type[n_plants:])
    assert n_cross == solution.objective
    # Crossings in the front
    assert n_cross < n_plants
    assert all(ty == "Node" for ty in solution.tree_type[:n_cross])
    assert all(ty == "Leaf" for ty in solution.tree_type[n_cross:])
    # Constraint on leaves
    assert all(
        solution.tree_left[i] == solution.tree_right[i] == 0
        for i in range(n_plants)
        if solution.tree_type[i] == "Leaf"
    )
    assert all(
        solution.tree_left[i] > 0
        and solution.tree_right[i] > 0
        and crossable(i)
        for i in range(n_plants)
        if solution.tree_type[i] == "Node"
    )
    # DAG
    # ideotype
    assert sum(is_ideotype(i) for i in range(n_plants)) == 1

    return all(
        [
            len(solution.tree_data) == n_tcells,
            len(solution.tree_type) == n_tcells,
            len(solution.tree_left) == n_tcells,
            len(solution.tree_right) == n_tcells,
            # Nulls at the back
            n_plants <= n_tcells,
            all(ty in ("Node", "Leaf") for ty in solution.tree_type[:n_plants]),
            all(ty == "Null" for ty in solution.tree_type[n_plants:]),
            n_cross == solution.objective,
            # Crossings in the front
            n_cross < n_plants,
            all(ty == "Node" for ty in solution.tree_type[:n_cross]),
            all(ty == "Leaf" for ty in solution.tree_type[n_cross:]),
            # Constraint on leaves
            all(
                solution.tree_left[i] == solution.tree_right[i] == 0
                for i in range(n_plants)
                if solution.tree_type[i] == "Leaf"
            ),
            all(
                solution.tree_left[i] > 0
                and solution.tree_right[i] > 0
                and crossable(i)
                for i in range(n_plants)
                if solution.tree_type[i] == "Node"
            ),
            # DAG
            # ideotype
            sum(is_ideotype(i) for i in range(n_plants)) == 1,
        ]
    )
