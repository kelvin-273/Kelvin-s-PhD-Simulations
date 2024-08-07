import eugene_rs
import eugene.utils as eu
from eugene.plant_models.plant2 import PlantSPC
from eugene.solution import BaseSolution
from typing import List, Optional


def breeding_program_distribute(
    dist_array: List[int], timeout=None
) -> Optional[BaseSolution]:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """
    n_loci = len(dist_array)
    pop_0 = eu.distribute_to_plants(dist_array)

    return breeding_program(n_loci, pop_0, timeout=timeout)


def breeding_program(
    n_loci: int, pop_0: List[PlantSPC], timeout=None
) -> Optional[BaseSolution]:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """

    result = eugene_rs.min_gen.naive1.breeding_program_python(
        n_loci,
        [
            [[bool(allele) for allele in row] for row in x.to_bitlist()]
            for x in pop_0
        ],
        timeout=timeout,
    )
    return BaseSolution(*result) if result else None


def mingen_answer(
    n_loci: int, pop_0: List[PlantSPC], timeout=None
) -> Optional[BaseSolution]:
    """
    Solves a distribute instance using the minizinc model.
    A MinizincContext can be passed in as an optional parameter,
    otherwise a MinizincContext is constructed using cp-sat from OR-Tools.
    """

    return eugene_rs.min_gen.naive1.mingen_answer(
        n_loci,
        [
            [[bool(allele) for allele in row] for row in x.to_bitlist()]
            for x in pop_0
        ],
        timeout=timeout,
    )
