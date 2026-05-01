from math import log, log2
from dataclasses import dataclass
from collections import Counter
from typing import NewType, Callable
from random import randrange

from eugene.utils import *

Chromosome = NewType("Chromosome", int)


@dataclass
class Plant:
    n_loci: int
    chrom1: int
    chrom2: int

    def format(self):
        xu = format(self.chrom1, f"0{self.n_loci}b")
        xl = format(self.chrom2, f"0{self.n_loci}b")
        return (xu, xl)

    def create_gamete(self) -> Chromosome:
        split = randrange(0, self.n_loci)
        start = randint(0, 1)

        return self.gamete_by_crosspoint(split, start)

    def gamete_by_crosspoint(self, start: int, split: int):
        """
        Creates a gamete with cross point
        """
        # bit-string representing the source of the ith allele ∀i∈[n]
        mask = (1 << self.n_loci + 1) - 1
        crossing = (mask >> split) << split
        if not start:
            crossing ^= mask

        return (crossing & self.chrom1) | (~crossing & self.chrom2)

    def create_gamete2(self, other: "Plant") -> Chromosome:
        split = randrange(0, self.n_loci)
        start = randint(0, 1)
        chroms = randrange(4)
        chrom_x_cross = self.chrom1 if chroms & 2 else self.chrom2
        chrom_x_clone = self.chrom2 if chroms & 2 else self.chrom1
        chrom_y_cross = other.chrom1 if chroms & 1 else other.chrom2
        chrom_y_clone = other.chrom2 if chroms & 1 else other.chrom1

        # bit-string representing the source of the ith allele ∀i∈[n]
        mask = (1 << self.n_loci + 1) - 1
        crossing = mask >> split << split
        if not start:
            crossing ^= mask

        return (
            chrom_x_clone,
            (crossing & chrom_x_cross) | (~crossing & chrom_y_cross),
            (~crossing & chrom_x_cross) | (crossing & chrom_y_cross),
            chrom_y_clone,
        )

    def cross(self, other):
        assert self.n_loci == other.n_loci
        gamete1 = self.create_gamete()
        gamete2 = other.create_gamete()

        return Plant(self.n_loci, gamete1, gamete2)

    def cross2(self, other: "Plant"):
        assert self.n_loci == other.n_loci
        gametes = self.create_gamete2(other)
        gamete1, gamete2 = sample(gametes, 2)

        return Plant(self.n_loci, gamete1, gamete2)

    def dom_strong(self, other):
        """
        Returns True iff every potential gamete of self dominates all
        potential gametes of other.
        """
        assert self.n_loci == other.n_loci
        return all(
            self.chrom1 | other.chrom1 == self.chrom1,
            self.chrom1 | other.chrom2 == self.chrom1,
            self.chrom2 | other.chrom1 == self.chrom2,
            self.chrom2 | other.chrom2 == self.chrom2,
        )

    def dom_weak(self, other: "Plant"):
        """
        Returns True iff the upper and lower chromosomes of other do not have
        favourable alleles that are not already in the respective chromosome of
        self.
        """
        assert self.n_loci == other.n_loci
        return (
            self.chrom1 | other.chrom1 == self.chrom1
            and self.chrom2 | other.chrom2 == self.chrom2
        )


    def reachable_gametes(self):
        """
        Returns the set of gametes that can be created from this plant.
        """
        d = Counter(
            self.gamete_by_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return set(d.keys())

    def reachable_gametes_with_counts(self):
        d = Counter(
            self.gamete_by_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return list(d.items())

    def __hash__(self):
        return (self.chrom1 << self.n_loci) | self.chrom2


@dataclass
class Population:
    plants: List[Plant]
    has_goal: bool


Strategy = Callable[[Population], Population]
ChoosingStrategy = Callable[[Population], Tuple[Plant, Plant]]


def generate_random_plant(n_loci: int) -> Plant:
    return Plant(n_loci, randrange(0, (1 << n_loci)), randrange(0, (1 << n_loci)))


def generate_goal(n_loci: int) -> Plant:
    g = (1 << n_loci) - 1
    return Plant(n_loci, g, g)


def union(pop: List[Plant]) -> Plant:
    assert len({x.n_loci for x in pop}) == 1
    acc = Plant(0, 0, 0)
    for x in pop:
        acc = Plant(x.n_loci, acc.chrom1 | x.chrom1, acc.chrom2 | x.chrom2)
    return acc


def prob_z_given_xy(z: Plant, x: Plant, y: Plant):
    assert z.n_loci == x.n_loci == y.n_loci
    nl = z.n_loci

    def aux(c: Chromosome, x: Plant):
        out = 0
        for start in [0, 1]:
            for split in range(nl):
                gx = x.gamete_by_crosspoint(start=start, split=split)
                if gx == c:
                    out += 1
        return out

    return aux(z.chrom1, x) * aux(z.chrom2, y) / (4 * nl**2)


def prob_z_given_xy_fast(z: Plant, x: Plant, y: Plant):
    assert z.n_loci == x.n_loci == y.n_loci
    nl = z.n_loci

    def len_matching_pre(s: str, t: str):
        out = 0
        for a, b in zip(s, t):
            if a != b:
                break
            out += 1
        return out

    def aux(c: Chromosome, x: Plant):
        c = format(c, f"0{nl}b")
        xu = format(x.chrom1, f"0{nl}b")
        xl = format(x.chrom2, f"0{nl}b")

        # length of prefix that matches
        lmpu = len_matching_pre(c, xu)
        lmpl = len_matching_pre(c, xl)
        lmsu = len_matching_pre(c[::-1], xu[::-1])
        lmsl = len_matching_pre(c[::-1], xl[::-1])
        matches_upper = max(0, min(nl - 1, lmpu) - (nl - lmsl) + 1)
        matches_lower = max(0, min(nl - 1, lmpl) - (nl - lmsu) + 1)
        return matches_upper + matches_lower

    return aux(z.chrom1, x) * aux(z.chrom2, y) / (4 * nl**2)


def number_of_trials_to_create(p, gamma):
    if p == 0:
        raise ValueError("genotype is unreachable")
    if p == 1:
        return 0
    return log(1 - gamma) / log(1 - p)
