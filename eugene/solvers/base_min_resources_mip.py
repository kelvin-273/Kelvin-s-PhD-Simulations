from typing import List
from eugene.solution import BaseSolution
from eugene.plant_models.plant2 import PlantSPC
from eugene.utils import distribute_to_plants
import gurobipy as gp
from gurobipy import GRB
from collections import defaultdict
from itertools import combinations
import eugene_pywrapper
from math import ceil


def build_model(pop, pop_0, target, H, crossings, c):
    m = gp.Model("min_res")

    pop = list(pop)
    pop_0 = set(pop_0)
    crossings = list(crossings)

    T = range(H)
    T_all = range(H + 1)

    incoming = defaultdict(list)
    for x, y, z in crossings:
        incoming[z].append((x, y, z))

    u = m.addVars(pop, T_all, vtype=GRB.BINARY, name="u")
    f = m.addVars(crossings, T, vtype=GRB.BINARY, name="f")

    m.setObjective(
        gp.quicksum(
            c[x, y, z] * f[x, y, z, t]
            for (x, y, z) in crossings
            for t in range(H)
        ),
        GRB.MINIMIZE,
    )

    m.addConstrs((u[z, 0] == 1 for z in pop_0), name="MIPMinResPopIn")
    m.addConstrs(
        (u[z, 0] == 0 for z in pop if z not in pop_0), name="MIPMinResPopOut"
    )
    m.addConstr(u[target, H] == 1, name="MIPMinResTarget")

    m.addConstrs(
        (u[x, t] >= f[x, y, z, t] for (x, y, z) in crossings for t in T),
        name="MIPMinResGameteRequirementL",
    )

    m.addConstrs(
        (u[y, t] >= f[x, y, z, t] for (x, y, z) in crossings for t in T),
        name="MIPMinResGameteRequirementR",
    )

    m.addConstrs(
        (
            u[z, t + 1]
            == u[z, t]
            + gp.quicksum(f[x, y, z, t] for (x, y, z) in incoming[z])
            for z in pop
            for t in T
        ),
        name="MIPMinResProduction",
    )

    return m, u, f


def build_model_2(pop, pop_0, target, H, crossings, c):
    m = gp.Model("min_res")

    pop = list(pop)
    pop_0 = set(pop_0)
    crossings = list(crossings)

    T = range(H)
    T_all = range(H + 1)

    incoming = defaultdict(list)
    for x, y, z in crossings:
        incoming[z].append((x, y, z))

    u = m.addVars(pop, T_all, vtype=GRB.BINARY, name="u")
    f = m.addVars(crossings, vtype=GRB.BINARY, name="f")

    m.setObjective(
        gp.quicksum(c[x, y, z] * f[x, y, z] for (x, y, z) in crossings),
        GRB.MINIMIZE,
    )

    m.addConstrs((u[z, 0] == 1 for z in pop_0), name="MIPMinResPopIn")
    m.addConstrs(
        (u[z, 0] == 0 for z in pop if z not in pop_0), name="MIPMinResPopOut"
    )
    m.addConstr(u[target, H] == 1, name="MIPMinResTarget")

    m.addConstrs(
        (u[z, t + 1] >= u[z, t] for z in pop for t in T),
        name="MIPMinResConservation",
    )

    m.addConstrs(
        (
            u[z, t + 1] + f[x, y, z] - u[x, t] <= 1
            for (x, y, z) in crossings
            for t in T
        ),
        name="MIPMinResGameteRequirementL",
    )

    m.addConstrs(
        (
            u[z, t + 1] + f[x, y, z] - u[y, t] <= 1
            for (x, y, z) in crossings
            for t in T
        ),
        name="MIPMinResGameteRequirementR",
    )

    m.addConstrs(
        (
            u[z, H]
            == u[z, 0] + gp.quicksum(f[x, y, z] for (x, y, z) in incoming[z])
            for z in pop
        ),
        name="MIPMinResProduction",
    )

    return m, u, f


def breeding_program_distribute(dist_array: List[int]) -> BaseSolution:
    n_loci = len(dist_array)
    pop_0 = distribute_to_plants(dist_array)
    return breeding_program(n_loci, pop_0)


class CrossingCost:
    def __init__(self, rec_rate: List[float], gamma: float):
        self.rec_rate = eugene_pywrapper.utils.PyRecRate(rec_rate)
        self.gamma = gamma

    def __getitem__(self, index):
        return self.rec_rate.cost_of_crossing(
            self.gamma, index[0], index[1], index[2]
        )


def breeding_program(
    n_loci: int,
    pop_0: List[PlantSPC],
    rec_rate: List[float],
    gamma: float,
    costs=None,
    timeout=None,
) -> BaseSolution:
    # Set a timeout for the whole function to prevent excessively long runtimes on large instances

    pop = list(range(1 << 2 * n_loci))
    pop_0 = {(x.chrom1 << n_loci) | x.chrom2 for x in pop_0}
    target = (1 << 2 * n_loci) - 1
    gen_bound = 2 * n_loci
    crossings = list(generate_crossings(n_loci))

    if costs is None:
        costs = CrossingCost(rec_rate, gamma)

    m, u, f = build_model(pop, pop_0, target, gen_bound, crossings, costs)

    # print("=== Built Model ===")
    m.optimize()
    if timeout is not None:
        m.setParam("TimeLimit", timeout)

    # print("=== Solution ===")
    # for x in pop:
    #     for t in range(gen_bound + 1):
    #         if u[x, t].X == 1:
    #             print(f"genotype {x} at gen {t}")

    # for x, y, z in crossings:
    #     for t in range(gen_bound):
    #         if f[x, y, z, t].X == 1:
    #             print(f"crossing ({x}, {y} -> {z}) at gen {t}")

    return BaseSolution(
        tree_data=[],
        tree_type=[],
        tree_left=[],
        tree_right=[],
        objective=ceil(m.ObjVal),
    )


def mincross_breeding_program(
    n_loci: int,
    pop_0: List[PlantSPC],
    rec_rate: List[float],
    gamma: float,
    costs=None,
    timeout=None,
):
    pop = list(range(1 << 2 * n_loci))
    pop_0 = {(x.chrom1 << n_loci) | x.chrom2 for x in pop_0}
    target = (1 << 2 * n_loci) - 1
    gen_bound = 2 * n_loci
    crossings = list(generate_crossings(n_loci))

    if costs is None:
        costs = CrossingCost(rec_rate, gamma)

    incoming = defaultdict(list)
    for x, y, z in crossings:
        incoming[z].append((x, y, z))

    m, u, f = build_model(pop, pop_0, target, gen_bound, crossings, costs)
    # m, u, f = build_model_2(pop, pop_0, target, gen_bound, crossings, costs)

    m.optimize()

    print("=== Built Model ===")
    m.optimize()
    if timeout is not None:
        m.setParam("TimeLimit", timeout)

    print("=== Solution ===")
    for x in pop:
        for t in range(gen_bound + 1):
            if u[x, t].X == 1:
                print(f"genotype {x} at gen {t}")

    for x, y, z in crossings:
        for t in range(gen_bound):
            if f[x, y, z, t].X == 1:
                print(f"crossing ({x}, {y} -> {z}) at gen {t}")

    genotype_idx_map = {}
    sol = BaseSolution(
        tree_data=[],
        tree_type=[],
        tree_left=[],
        tree_right=[],
        objective=ceil(m.ObjVal),
    )

    def dfs(z, t):
        if z in genotype_idx_map:
            return
        if t == 0:
            # z is a leaf
            assert u[z, t].X == 1
            genotype_idx_map[z] = len(sol.tree_data)
            sol.tree_data.append(genotype_int_to_lists(z, n_loci))
            sol.tree_type.append("Leaf")
            sol.tree_left.append(0)
            sol.tree_right.append(0)
        elif u[z, t].X == u[z, t-1].X:
            dfs(z, t-1)
        else:
            assert u[z, t].X == 1 and u[z, t-1].X == 0
            # find crossing
            for x, y in generate_crossings_for_offpring(n_loci, z):
                if f[x, y, z, t-1].X == 1:
                    dfs(x, t-1)
                    dfs(y, t-1)
                    genotype_idx_map[z] = len(sol.tree_data)
                    sol.tree_data.append(genotype_int_to_lists(z, n_loci))
                    sol.tree_type.append("Node")
                    sol.tree_left.append(genotype_idx_map[x])
                    sol.tree_right.append(genotype_idx_map[y])
                    return
            else:
                raise ValueError(f"Could not find crossing for {z} at gen {t}")

    dfs(target, gen_bound)

    # Reverse the tree to ensure that parents come before children
    n_cells = len(sol.tree_data)
    for i in range(n_cells // 2):
        # Swap i and n_cells - 1 - i
        sol.tree_data[i], sol.tree_data[n_cells - 1 - i] = sol.tree_data[n_cells - 1 - i], sol.tree_data[i]
        sol.tree_type[i], sol.tree_type[n_cells - 1 - i] = sol.tree_type[n_cells - 1 - i], sol.tree_type[i]
        sol.tree_left[i], sol.tree_left[n_cells - 1 - i] = sol.tree_left[n_cells - 1 - i], sol.tree_left[i]
        sol.tree_right[i], sol.tree_right[n_cells - 1 - i] = sol.tree_right[n_cells - 1 - i], sol.tree_right[i]
    for i in range(n_cells):
        if sol.tree_type[i] == "Node":
            # Update left and right indices to reflect the reversal
            sol.tree_left[i] = n_cells - 1 - sol.tree_left[i]
            sol.tree_right[i] = n_cells - 1 - sol.tree_right[i]

    return sol


def genotype_int_to_lists(z: int, n_loci: int):
    z_list = [[], []]
    mask = (1 << (2 * n_loci)) - 1
    z_1 = (z >> n_loci) & mask
    z_2 = z & mask
    for i in range(n_loci):
        z_list[0].append((z_1 >> i) & 1)
        z_list[1].append((z_2 >> i) & 1)
    return z_list


def mincross_breeding_program_distribute(
    dist_array: List[int],
    rec_rate: List[float],
    gamma: float,
    timeout=None,
    costs=None,
):
    n_loci = len(dist_array)
    pop_0 = distribute_to_plants(dist_array)
    if costs is None:
        costs = CrossingCost(rec_rate, gamma)
    return mincross_breeding_program(n_loci, pop_0, rec_rate, gamma, costs)


def upper_bound(n_loci: int, pop_0: List[PlantSPC]) -> int:
    return n_loci


def generate_crossings(n_loci):
    for x, z_1 in generate_recombinations(n_loci):
        for y, z_2 in generate_recombinations(n_loci):
            z = (z_1 << n_loci) | z_2
            yield x, y, z


def generate_crossings_for_offpring(n_loci, z):
    mask = (1 << n_loci) - 1
    for x_1, x_2 in generate_masks(n_loci):
        z_1 = z >> n_loci
        x_1 = mask ^ z_1 ^ x_1
        x_2 = mask ^ z_1 ^ x_2
        for y_1, y_2 in generate_masks(n_loci):
            z_2 = z & mask
            y_1 = mask ^ z_2 ^ y_1
            y_2 = mask ^ z_2 ^ y_2
            x = (x_1 << n_loci) | x_2
            y = (y_1 << n_loci) | y_2
            yield x, y


def generate_recombinations(n_loci: int):
    mask = (1 << n_loci) - 1
    for z_1 in range(1 << n_loci):
        for x_1_mask, x_2_mask in generate_masks(n_loci):
            x_1 = mask ^ z_1 ^ x_1_mask
            x_2 = mask ^ z_1 ^ x_2_mask
            x = (x_1 << n_loci) | x_2
            yield x, z_1


def generate_masks(n_loci: int):
    mask = (1 << n_loci) - 1
    yield (mask, mask)
    for m in range(1, n_loci + 1):
        for gaps in combinations(range(n_loci), r=m):
            z = mask
            for j in gaps:
                z ^= 1 << j
            for k in gaps:
                mask_r = mask >> (n_loci - k)
                mask_l = mask ^ mask_r
                z_1 = mask_l | z
                z_2 = mask_r | z
                yield z_1, z_2
                yield z_2, z_1
