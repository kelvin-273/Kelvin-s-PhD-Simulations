from collections import Counter, defaultdict
from dataclasses import dataclass
from random import randint, randrange, random, sample
from typing import NewType, List
from abc import ABC, abstractmethod, abstractclassmethod
from bitarray import frozenbitarray
from bitarray.util import urandom, zeros


Chromosome = NewType("Chromosome", int)


@dataclass
class Crossable(ABC):
    """Abstract class for plants"""

    @abstractmethod
    def crosspoints(self):
        pass

    @abstractmethod
    def random_crosspoint(self):
        pass

    @abstractmethod
    def gamete_specified(self, crosspoint):
        pass

    def gamete_random(self):
        crosspoint = self.random_crosspoint()
        return self.gamete_specified(crosspoint)

    @abstractmethod
    def cross_specified(self, other, crosspoints):
        pass

    def cross_random(self, other):
        self.cross_specified(
            other, (self.random_crosspoint(), other.random_crosspoint())
        )

    @abstractmethod
    def from_gametes(*args):
        pass


class DomWeak(ABC):
    """Abstract class with weak dominance relation"""

    @abstractmethod
    def dom_weak(self, other) -> bool:
        pass


class DomStrong(DomWeak):
    """Abstract class with strong dominance relation"""

    @abstractmethod
    def dom_strong(self, other) -> bool:
        pass


class Unionable(ABC):
    """Abstract class for feasibility checks"""

    @abstractclassmethod
    def union(pop):
        from sys import stderr

        stderr.write(
            "DeprecationWarning: union will be replaced with is_feasible\n"
        )
        # raise DeprecationWarning("union will be replaced with is_feasible")

    @abstractclassmethod
    def is_feasible(pop) -> bool:
        pass

    @abstractclassmethod
    def goal(n_loci):
        pass


class InitialPop(ABC):
    """Abstract class for initial population generators"""

    @abstractmethod
    def initial_pop_random(n_loci, n_elements):
        pass

    @abstractmethod
    def initial_pop_trait_introgression(n_loci, n_holes, n_elite, n_donor):
        pass

    @abstractmethod
    def initial_pop_singles_homo(n_loci):
        pass

    @abstractmethod
    def initial_pop_singles_hetero(n_loci):
        pass


@dataclass(order=True)
class PlantSPC(Crossable, DomStrong):
    """Plant with single-point crossover"""

    n_loci: int
    chrom1: int
    chrom2: int

    def from_gametes(n_loci, chrom1, chrom2):
        return PlantSPC(n_loci, chrom1, chrom2)

    def cross_specified(self, other, crosspoints):
        k1, k2 = crosspoints
        g1 = self.gamete_specified(k1)
        g2 = other.gamete_specified(k2)
        return PlantSPC(self.n_loci, g1, g2)

    def random_crosspoint(self):
        """
        Generates a random crosspoint as a tuple (start, split) where
        start in {0, 1} is the chromosome whose prefix forms the prefix of the
        gamete and split in {0..n-1} (inclusive) is the locus up to which the prefix of
        the start chromosome is created from.
        Mainly for use with random gamete creation.
        """
        return (randint(0, 1), randrange(0, self.n_loci))

    def gamete_specified(self, crosspoint):
        """
        Creates a gamete with crosspoint = (start, split)
        """
        assert type(crosspoint) is tuple and len(crosspoint) == 2
        assert type(crosspoint[0]) is int and type(crosspoint[1]) is int
        start, split = crosspoint
        # bit-string representing the source of the ith allele ∀i∈[n]
        mask = (1 << self.n_loci + 1) - 1
        crossing = (mask >> split) << split
        if not start:
            crossing ^= mask

        return (crossing & self.chrom1) | (~crossing & self.chrom2)

    def cross_random(self, other):
        return PlantSPC(
            self.n_loci, self.gamete_random(), other.gamete_random()
        )

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

    def dom_weak(self, other):
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
            self.create_gamete_with_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return set(d.keys())

    def reachable_gametes_with_counts(self):
        d = Counter(
            self.create_gamete_with_crosspoint(i, j)
            for i in range(2)
            for j in range(self.n_loci)
        )
        return list(d.items())

    def crosspoints(self):
        for i in range(2):
            for j in range(self.n_loci):
                yield (i, j)

    def __hash__(self):
        return (self.chrom1 << self.n_loci) | self.chrom2

    @staticmethod
    def initial_pop_trait_introgression(
        n_loci: int, n_holes: int, n_elite=1, n_donor=1
    ):
        """
        Creates an initial population for trait introgression.
        Output list is divided into [elite_pop : donor_pop].
        Elite and donor populations each have a complimentary set of homozygous
        loci and loci not in that set are hetrozygous or homozygous
        unfavourable.
        """
        assert n_holes <= n_loci // 2
        holes = sample(range(n_loci), n_holes)

        bitbed = (1 << n_loci) - 1
        elite_mask = bitbed
        for i in holes:
            elite_mask ^= 1 << i

        def gen_elite(mask_elite: int):
            xu = xl = (1 << n_loci) - 1
            # get bits that need fixing
            d = (~mask_elite) & xu & xl
            while d:
                u = randrange(1 << n_loci)
                l = randrange(1 << n_loci)
                xu = ((~d) & xu) | (d & u)
                xl = ((~d) & xl) | (d & l)
                d = (~mask_elite) & xu & xl
            return PlantSPC(n_loci, xu, xl)

        elite_pop = [gen_elite(elite_mask) for _ in range(n_elite)]
        donor_pop = [gen_elite(bitbed ^ elite_mask) for _ in range(n_elite)]

        return elite_pop + donor_pop

    def initial_pop_trait_introgression_homo(
        n_loci: int, n_holes: int, n_elite=1, n_donor=1
    ):
        """
        Creates an initial population for trait introgression.
        Output list is divided into [elite_pop : donor_pop].
        Elite and donor populations each have a complimentary set of homozygous
        loci and loci not in that set are hetrozygous or homozygous
        unfavourable.
        """
        assert n_holes <= n_loci // 2
        holes = sample(range(n_loci), n_holes)

        bitbed = (1 << n_loci) - 1
        elite_mask = bitbed
        for i in holes:
            elite_mask ^= 1 << i
        donor_mask = bitbed ^ elite_mask

        elite_pop = [PlantSPC(n_loci, elite_mask, elite_mask)]
        donor_pop = [PlantSPC(n_loci, donor_mask, donor_mask)]

        return elite_pop + donor_pop

    @staticmethod
    def initial_pop_random(n_loci: int, n_individuals: int, p=0.5):
        """
        Generates a population of uniform randomly generated plants such that
        the ideotype is in its span.
        """
        mask_max = (1 << n_loci) - 1

        def gen_chrom(n_loci, p):
            c = 0
            for _ in range(n_loci):
                c <<= 1
                c |= random() < p
            return c

        def union(pop):
            out = 0
            for a in pop:
                out |= a
            return out

        pop = [gen_chrom(n_loci, p) for _ in range(2 * n_individuals)]
        mask_cur = union(pop)
        while mask_cur != mask_max:
            pop = [
                x | (gen_chrom(n_loci, p) & (mask_max ^ mask_cur)) for x in pop
            ]
            mask_cur = union(pop)

        return [
            PlantSPC(n_loci, pop[2 * i], pop[2 * i + 1])
            for i in range(n_individuals)
        ]

    @staticmethod
    def initial_pop_singles_homo(n_loci: int):
        return [PlantSPC(n_loci, 1 << i, 1 << i) for i in range(n_loci)]

    @staticmethod
    def initial_pop_singles_hetero(n_loci: int):
        return [PlantSPC(n_loci, 1 << i, 0) for i in range(n_loci)] + [
            PlantSPC(n_loci, 0, 1 << i) for i in range(n_loci)
        ]

    @staticmethod
    def initial_pop_distribute(xs: list[int]):
        """
        Generates the initial population for a distribute instance.
        """
        n_loci = len(xs)
        pop = defaultdict(lambda: PlantSPC(n_loci, 0, 0))
        for i, x in enumerate(xs):
            x = pop[i]
            x.chrom1 |= 1 << (n_loci - i - 1)
            x.chrom2 |= 1 << (n_loci - i - 1)
        return list(pop.values())

    @staticmethod
    def union(pop) -> int:
        """
        Takes the union of all chromosomes in the population.
        Useful for testing feasibility.

        :pop: List of individual PlantSPC's
        :returns: integer containing the bits after union of all chromosomes in
            the population
        """
        out = 0
        for x in pop:
            out |= x.chrom1 | x.chrom2
        return out

    def __str__(self):
        s1 = format(self.chrom1, f"0{self.n_loci}b")
        s2 = format(self.chrom2, f"0{self.n_loci}b")
        return f"({s1}, {s2})"

    def to_bitlist(self):
        upper = self.chrom1
        lower = self.chrom2

        out = [[0] * self.n_loci, [0] * self.n_loci]
        for i in reversed(range(self.n_loci)):
            upper, out[0][i] = upper >> 1, upper & 1
            lower, out[1][i] = lower >> 1, lower & 1

        return out


@dataclass(order=True)
class PlantSPCBitarray:
    """Plant with single-point crossover"""

    n_loci: int
    chrom1: frozenbitarray
    chrom2: frozenbitarray

    @staticmethod
    def initial_pop_random(n_loci: int, n_individuals: int, p=0.5):
        """
        Generates a population of uniform randomly generated plants such that
        the ideotype is in its span.
        """
        mask_max = ~zeros(n_loci)

        def gen_chrom(n_loci, p):
            c = zeros(n_loci)
            for i in range(n_loci):
                c[i] = random() < p
            return c

        def union(pop):
            out = zeros(n_loci)
            for a in pop:
                out |= a
            return out

        pop = [gen_chrom(n_loci, p) for _ in range(2 * n_individuals)]
        mask_cur = union(pop)

        while mask_cur != mask_max:
            pop = [x | (gen_chrom(n_loci, p) & ~mask_cur) for x in pop]
            mask_cur = union(pop)
        return [
            PlantSPCBitarray(n_loci, pop[2 * i], pop[2 * i + 1])
            for i in range(n_individuals)
        ]

    @staticmethod
    def initial_pop_distribute(xs: list[int]):
        """
        Generates the initial population for a distribute instance.
        """
        n_loci = len(xs)
        pop = defaultdict(lambda: PlantSPCBitarray(n_loci, zeros(n_loci), zeros(n_loci)))
        for i, x in enumerate(xs):
            x = pop[i]
            x.chrom1[i] = 1
            x.chrom2[i] = 1
        return list(pop.values())

    @staticmethod
    def union(n_loci, pop) -> int:
        """
        Takes the union of all chromosomes in the population.
        Useful for testing feasibility.

        :pop: List of individual PlantSPC's
        :returns: integer containing the bits after union of all chromosomes in
            the population
        """
        out = zeros(n_loci)
        for x in pop:
            out |= x.chrom1 | x.chrom2
        return frozenbitarray(out)

    def to_bitlist(self):
        upper = self.chrom1
        lower = self.chrom2
        return [upper.tolist(), lower.tolist()]

    def to_ints(self):
        upper = self.chrom1
        lower = self.chrom2
        return int(upper.to01(), 2), int(lower.to01(), 2)

    def to_PlantSPC(self):
        chrom1, chrom2 = self.to_ints()
        return PlantSPC(self.n_loci, chrom1, chrom2)

    def to_bool_list(self):
        upper = self.chrom1
        lower = self.chrom2
        return [
            [b == 1 for b in upper.tolist()],
            [b == 1 for b in lower.tolist()],
        ]


@dataclass
class WDataG(DomStrong, Unionable):
    """Docstring for WData."""

    def __init__(self, x, *, history):
        self.x = x
        self.history = history

    def dom_gamete(self, other):
        return not (other.x & ~self.x)

    def dom_weak(self, other):
        return self.x.dom_weak(other.x)

    def dom_strong(self, other):
        return self.x.dom_weak(other.x)

    def union(pop):
        raise Exception("not implemented")

    def is_feasible(pop):
        raise Exception("not implemented")

    def goal(n_loci):
        return (1 << n_loci) - 1

    def __hash__(self):
        return self.x.__hash__()

    def __invert__(self):
        return WDataG(~self.x, history=self.history)


@dataclass
class WDataP(Crossable, DomStrong, Unionable):
    """Docstring for WData."""

    def __init__(self, x, *, history, count=1):
        self.x = x
        self.history = history

    @property
    def n_loci(self):
        return self.x.n_loci

    # TODO: Dependence on fixed ploidy level should be removed <05-07-22> #
    @property
    def chrom1(self):
        return self.x.chrom1

    @property
    def chrom2(self):
        return self.x.chrom2

    def crosspoints(self):
        return self.x.crosspoints()

    def random_crosspoint(self):
        return self.x.random_crosspoint()

    def gamete_specified(self, crosspoint):
        return WDataG(
            self.x.gamete_specified(crosspoint),
            history=self.history,
            count=self.count,
        )

    def gamete_random(self):
        return WDataG(
            self.x.gamete_random(),
            history=self.history,
            count=self.count,
        )

    def cross_specified(self, other, crosspoints):
        return WDataP(
            self.x.cross_specified(other.x, crosspoints),
            history=(self.history, other.history),
        )

    def cross_random(self, other):
        return self.cross_specified(self, other, self.random_crosspoint())

    def from_gametes(n_loci, chrom1, chrom2):
        return WDataP(
            PlantSPC(n_loci, chrom1.x, chrom2.x),
            history=(chrom1.history, chrom2.history),
        )

    def dom_weak(self, other):
        return self.x.dom_weak(other.x)

    def dom_strong(self, other):
        return self.x.dom_strong(other.x)

    @classmethod
    def union(pop):
        raise Exception("not implemented")

    @classmethod
    def is_feasible(pop):
        raise Exception("not implemented")

    def goal(n_loci):
        raise DeprecationWarning(
            """this might not be necessary given that we set the goal in the
            breeding program class"""
        )
        return WDataP(
            PlantSPC(n_loci, (1 << n_loci) - 1, (1 << n_loci) - 1),
            history=None,
        )

    # TODO: Is there a way to automatically derive this <05-07-22> #
    def __hash__(self):
        return hash((self.chrom1 << self.n_loci) & self.chrom2)


def wd_plant(plant_model):
    class WD(plant_model):
        """Class that contains the historical data of the plant"""

        def __init__(self, *args, **kwargs):
            """TODO: to be defined."""
            plant_model.__init__(self, *args, **kwargs)
            self.histories = None

    return WD


if __name__ == "__main__":
    x = PlantSPC(4, 5, 10)
    wx = WDataP(x, history=1)
    print(wx.random_crosspoint())
    print(wx)
    print(wx.x)
    print(wx.n_loci)
    print(wx.chrom1)
    print(wx.chrom2)
