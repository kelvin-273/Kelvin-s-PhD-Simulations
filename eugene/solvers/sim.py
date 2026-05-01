"""
Plants
- Genomes
"""

from typing import Optional, Callable
from dataclasses import dataclass
from random import randrange, sample, choices
from collections import namedtuple

from eugene.plant_models.plant import (
    Plant,
    Population,
    generate_random_plant,
    generate_goal,
)
from eugene.plant_models.plant import (
    union,
    prob_z_given_xy_fast,
    number_of_trials_to_create,
)


@dataclass
class Args:
    n_loci: int
    n_remaining_loci: Optional[int]
    n_initial_pop: int
    gamma: float
    choose_parents: Callable
    choose_intermediate_target: Callable
    generate_initial_pop: Callable
    resource_cap: Optional[int] = None
    generation_cap: Optional[int] = None
    pruning: bool = False


def seq_breeding(args: Args):
    n_loci = args.n_loci
    choose_parents = args.choose_parents
    choose_intermediate_target = args.choose_intermediate_target
    generate_initial_pop = args.generate_initial_pop
    resource_cap = args.resource_cap
    generation_cap = args.generation_cap
    pruning = args.pruning
    gamma = args.gamma

    Results = namedtuple(
        "Results", ["n_generations", "n_plants_max", "n_plants_tot", "success"]
    )

    pop_0 = generate_initial_pop(args)
    if not union(pop_0) == generate_goal(n_loci):
        return Results(0, 0, 0, False)

    pop = pop_0.copy()
    goal = generate_goal(n_loci)
    t = 0
    n_tot = 0
    n_max = 0

    if pruning:
        pop = filter_non_dominating(pop)

    while goal not in pop and (not generation_cap or t < generation_cap):
        x, y = choose_parents(pop)
        z = choose_intermediate_target(x, y)

        if pruning and any(x0.dom_weak(z) for x0 in pop):
            continue
        # remove plants in the population that are dominated by z
        pop = [x for x in pop if not z.dom_weak(x)] + [z]

        p = prob_z_given_xy_fast(z, x, y)
        n = number_of_trials_to_create(p, gamma)

        if resource_cap and n > resource_cap:
            continue

        # all of the updating
        n_max = max(n, n_max)
        n_tot += n
        t += 1

    return Results(
        n_generations=t,
        n_plants_max=n_max,
        n_plants_tot=n_tot,
        success=goal in pop,
    )


def pcv_breeding(args: Args):
    pass


def filter_non_dominating(pop: Population) -> Population:
    n = len(pop)
    to_keep = [True] * n
    for i, x in enumerate(pop):
        for j in range(i):
            if not to_keep[j]:
                continue
            y = pop[j]
            if y.dom_weak(x):
                to_keep[i] = False
                break
            if x.dom_weak(y):
                to_keep[j] = False
    return [pop[i] for i in range(n) if to_keep[i]]


class PopulationGenerators:
    @staticmethod
    def generate_initial_pop_trait_introgression(args: Args):
        n_initial_pop = args.n_initial_pop
        n_loci = args.n_loci
        n_remaining_loci = args.n_remaining_loci
        assert n_remaining_loci <= n_loci // 2
        holes = sample(range(n_loci), n_remaining_loci)

        mask = (1 << n_loci) - 1
        for i in holes:
            mask ^= 1 << i

        def gen_elite(mask: int):
            xu = xl = (1 << n_loci) - 1
            d = (~mask) & xu & xl
            while d:
                u = randrange(1 << n_loci)
                l = randrange(1 << n_loci)
                xu = ((~d) & xu) | (d & u)
                xl = ((~d) & xl) | (d & l)
                d = (~mask) & xu & xl
            return Plant(n_loci, xu, xl)

        elite_pop = [
            # Plant(n_loci, randrange(1 << n_loci) | mask, randrange(1 << n_loci) | mask,)
            gen_elite(mask)
            for _ in range(n_initial_pop - 1)
        ]

        donor_pop = [gen_elite(~mask) for _ in range(1)]
        return elite_pop + donor_pop

    @staticmethod
    def generate_initial_pop_random(args: Args):
        n_initial_pop = args.n_initial_pop
        n_loci = args.n_loci
        return [generate_random_plant(n_loci) for _ in range(n_initial_pop)]

    @staticmethod
    def generate_initial_pop_single(args: Args):
        n_loci = args.n_loci
        n_remaining_loci = args.n_remaining_loci
        holes = sample(range(n_loci), n_remaining_loci)
        mask = 1 << n_loci
        mask -= 1
        chrom1 = chrom2 = mask
        for i in holes:
            if randrange(2):
                chrom1 ^= 1 << i
            else:
                chrom2 ^= 1 << i
        return [Plant(n_loci, chrom1, chrom2)]


class SelectionMethods:
    @staticmethod
    def choose_parents_by_ranking_method(f: Callable) -> Callable:
        pass

    @staticmethod
    def choose_parents_proportionate_to_score(f: Callable) -> Callable:
        """TODO: Docstring for choose_parents_proportionate_to_score.

        :f: TODO
        :returns: TODO

        """
        pass

    @staticmethod
    def GEBV(x: Plant) -> int:
        """
        Genomic Estimated Breeding Value under uniform weighted loci
        TODO: with Python3.10 this is made way simpler with .bit_count()
        """
        out = 0
        z = x.chrom1
        while z > 0:
            z, r = z >> 1, z & 1
            out += r

        z = x.chrom2
        while z > 0:
            z, r = z >> 1, z & 1
            out += r

        return out

    @staticmethod
    def OHV(x: Plant) -> int:
        """
        Optimal Haploid Value
        """
        out = 0
        c1, c2 = x.chrom1, x.chrom2
        while c1 > 0 or c2 > 0:
            c1, r1 = c1 >> 1, c1 & 1
            c2, r2 = c2 >> 1, c2 & 1
            out += r1 | r2
        return out

    @staticmethod
    def PCV(x: Plant, y: Plant) -> float:
        """
        Predicted Cross Value
        computes the probability that a gamete produced by a child of x and y
        will consist of only favourable alleles.

        In Han et al. 2017, this requires a recombination vector to determine
        the probability of generating a particular gamete because gamete
        production was performed using multipoint crossover.
        """
        pass

    @staticmethod
    def PCV_approx(x: Plant, y: Plant, n_trials=1000) -> float:
        assert x.n_loci == y.n_loci
        all_ones = (1 << x.n_loci) - 1
        count = 0  # frequentist
        for _ in range(n_trials):
            z = x.cross(y)
            g = z.create_gamete()
            count += g == all_ones
        return count / n_trials

    @staticmethod
    def choose_parents_uniform(pop: Population):
        return sample(pop, 2)

    @staticmethod
    def choose_parents_OHV(pop: Population):
        n = len(pop)
        p1 = max(range(n), key=lambda i: SelectionMethods.OHV(pop[i]))
        temp_pop = pop[:p1] + pop[p1 + 1:]
        p2 = max(range(n - 1), key=lambda i: SelectionMethods.OHV(temp_pop[i]))
        if p2 >= p1:
            p2 += 1
        return pop[p1], pop[p2]

    @staticmethod
    def choose_parents_GEBV(pop: Population):
        n = len(pop)
        p1 = max(range(n), key=lambda i: SelectionMethods.GEBV(pop[i]))
        temp_pop = pop[:p1] + pop[p1 + 1:]
        p2 = max(
            range(n - 1), key=lambda i: SelectionMethods.GEBV(temp_pop[i])
        )
        if p2 >= p1:
            p2 += 1
        return pop[p1], pop[p2]


def choose_parents_uniform(pop: Population):
    return choices(pop, k=2)


def choose_intermediate_target_uniform(x: Plant, y: Plant) -> Plant:
    return x.cross(y)


def simulation(args: Args):
    import matplotlib.pyplot as plt
    import seaborn as sn
    import numpy as np

    n_loci = args.n_loci
    gamma = args.gamma
    n_remaining_loci = args.n_remaining_loci
    xs = [seq_breeding(args) for _ in range(1000)]
    sn.scatterplot(
        [res.n_generations for res in xs], [res.n_plants_max for res in xs]
    )
    plt.title(
        "max_plants per generation vs generations to find target"
        + f"\nn_loci: {n_loci}"
        + f"\n$\gamma$: {gamma}"
        + f"\nelite alleles: {n_loci - n_remaining_loci}, donor alleles: {n_remaining_loci}"
    )
    plt.xlabel("no. generations")
    plt.ylabel("max. progeny per gen")
    plt.show()
    print("mean:", np.mean([res.n_generations for res in xs]))
    print("std:", np.std([res.n_generations for res in xs]))
    print("%success:", np.mean([res.success for res in xs]))


if __name__ == "__main__":
    x = Plant(5, 17, 24)
    y = Plant(5, 7, 30)
