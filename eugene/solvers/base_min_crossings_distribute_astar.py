import sys
import dataclasses
import heapq
import eugene.solvers.base_min_crossings_minizinc as em
from typing import List, Optional, Tuple
from eugene.solution import BaseSolution
from eugene_rs import min_cross


def breeding_program_distribute(
    dist_array: List[int], timeout=None, debug=False
) -> Optional[BaseSolution]:
    res = min_cross.distribute_astar.breeding_program_distribute_general_python(
        dist_array, timeout=timeout, diving=False, full_join=True,
        dominance=False
    )
    # res = astar(dist_array, ctx=BRANCHING_CTX, debug=debug)
    obj = res["objective"]
    # obj = ceil(res.f / 2)
    return BaseSolution(
        tree_data=[], tree_type=[], tree_left=[], tree_right=[], objective=obj
    )


BRANCHING_CTX = em.MinizincContext.from_solver_and_model_file(
    "gecode", "./eugene/solvers/minizinc/distribute_recomb.mzn"
)


@dataclasses.dataclass
class Node:
    xs: List[int]
    parent_gametes: Optional[Tuple[int, int]]
    parent_node: Optional[object]
    n_gametes: int
    n_segments: int
    g: int
    f: int

    def __lt__(self, other):
        return self.f < other.f

    @staticmethod
    def from_dist_array(dist_array):
        n_gametes = max(dist_array) + 1
        n_segments = len(dist_array)
        return Node(
            xs=dist_array,
            parent_gametes=None,
            parent_node=None,
            n_gametes=n_gametes,
            n_segments=n_segments,
            g=0,
            f=n_segments + n_gametes,
        )

    def __str__(self):
        return f"{self.xs} -> ({self.n_gametes}, {self.n_segments})"

    def _repr__(self):
        return (
            "<"
            + ", ".join(
                [
                    f"xs: {self.xs}",
                    f"parent: {None if self.parent_node is None else self.parent_node.xs}",
                    f"g: {self.g}, f: {self.f}",
                ]
            )
            + ">"
        )


def print_node(node: Node, file=None):
    print(
        "- {{{}}}".format(
            ", ".join(
                [
                    "type: node",
                    f"id: '{node.xs}'",
                    f"pId: '{node.parent_node.xs if node.parent_node else 'null'}'",
                    f"f: {node.f}",
                    f"g: {node.g}",
                    f"h: {node.n_gametes + node.n_segments}",
                ]
            )
        ),
        file=file if file else sys.stdout,
    )
    pass


def print_expansion(zs, xs, file=None):
    print(
        "- {{{}}}".format(
            ", ".join(
                [
                    "type: successor",
                    f"id: '{zs}'",
                    f"pId: '{xs}'",
                ]
            )
        ),
        file=file if file else sys.stdout,
    )
    pass


def astar(dist_array, ctx, debug=False):
    open_list = [Node.from_dist_array(dist_array)]
    open_dict = {str(open_list[0]): open_list[0]}
    closed_list = set()
    print("Decrease key needs to be implemented", file=sys.stderr)

    f = open(f"posthoc-{''.join(map(str, dist_array))}.trace.yaml", "w")
    print("version: 1.4.0", file=f)
    print("events:", file=f)

    while len(open_list) > 0:

        node = heapq.heappop(open_list)
        # print(f"- {{ type: node, id: \"{node.xs}\", pId: \"{node.parent_node.xs if node.parent_node else 'null'}\"}}", file=f)
        print_node(node, file=f)

        if debug:
            print(node.xs)

        if str(node) in closed_list:
            continue
        closed_list.add(str(node))
        open_dict.pop(str(node))

        if success(node):
            return node

        # TODO: is there something about updating costs that I'm missing

        children = branching(node, ctx=ctx)
        for child in children:
            key = str(child)
            print_expansion(child.xs, node.xs, file=f)
            if key not in closed_list and key not in open_dict:
                heapq.heappush(open_list, child)
                open_dict[key] = child
    f.close()
    return None


def branching(state: Node, ctx) -> List[Node]:
    first_full_join_solution = first_full_join_raw(state.xs)
    # first_full_join_solution = None
    if first_full_join_solution is not None:
        zs, gx, gy = first_full_join_solution
        zs = _simplify_dist_array(zs, len(zs) - 1)
        n_segments = len(zs)
        n_gametes = len(set(zs))
        return [
            Node(
                xs=zs,
                parent_gametes=(gx, gy),
                parent_node=state,
                n_gametes=n_gametes,
                n_segments=n_segments,
                g=state.g + 2,
                f=state.g + 2 + n_gametes + n_segments,
            )
        ]
    zss = generate_simplified_redistributions_brute_force(state.xs)
    return [
        Node(
            xs=zs,
            parent_gametes=(0, 1),
            parent_node=state,
            n_gametes=max(zs) + 1,
            n_segments=len(zs),
            g=state.g + 2,
            f=state.g + 2 + len(zs) + max(zs) + 1,
        )
        for zs in zss
    ]
    with ctx.instance.branch() as instance:
        instance["instance"] = state.xs
        instance["nLoci"] = state.n_segments
        instance["nPop"] = state.n_gametes
        all_solutions = instance.solve(all_solutions=True)

        # filter for unique solutions
        d = {}
        for res in all_solutions:
            print(f"res: {res}")
            zs = _simplify_dist_array(res.xs, res.nSegments)
            print(f"zs: {zs}")
            key = str(zs)
            if key not in d:
                d[key] = Node(
                    xs=_simplify_dist_array(res.xs, res.nSegments),
                    parent_gametes=(res.gx, res.gy),
                    parent_node=state,
                    n_gametes=res.nGametes,
                    n_segments=res.nSegments,
                    g=state.g + 2,
                    f=state.g + 2 + res.nGametes + res.nSegments,
                )

        # print(all_solutions)
        return list(d.values())


def success(state: Node):
    return len(state.xs) > 0 and all(x == state.xs[0] for x in state.xs)


def _simplify_dist_array(xs, n_segments):
    xs = xs.copy()  # get rid of this for more performance (and instability)
    n_loci = len(xs)
    mapping = [None] * (max(xs) + 1)  # assuming that xs is 0-indexed
    prev_element = xs[0]
    mapping[xs[0]] = 0
    xs[0] = mapping[xs[0]]
    x_max = 0
    i = j = 1
    while j < n_loci:
        current_element = xs[j]
        if xs[j] != prev_element:
            if mapping[xs[j]] is not None:
                xs[i] = mapping[xs[j]]
            else:
                x_max += 1
                mapping[xs[j]] = x_max
                xs[i] = x_max
            i += 1
        prev_element = current_element
        j += 1
        # print(f"xs: {xs},\tmapping: {mapping},\ti: {i},\tj: {j}")
    return xs[:i]


def first_full_join_raw(xs):
    n_loci = len(xs)
    n_pop = max(xs) + 1  # assumes xs is 0-indexed

    start_points = [n_loci - 1] * n_pop
    end_points = [0] * n_pop
    for i, gx in enumerate(xs):
        end_points[gx] = i
    for i, gx in reversed(list(enumerate(xs))):
        start_points[gx] = i

    # find the (gx,gy) with the smallest start_point[gx] such that
    # a full join occurs between gx and gy
    for gy, sy in enumerate(start_points):
        if sy > 0 and end_points[xs[sy - 1]] == sy - 1:
            gx = xs[sy - 1]
            out = xs.copy()
            for i in range(sy, end_points[gy] + 1):
                if xs[i] == gy:
                    out[i] = gx
            return (out, gx, gy)
    return None


def is_subsequence_of(xs, ys):
    return False


def is_distribute(xs):
    """
    Excludes the "no adjacent loci are owned by the same gamete" rule.
    """
    if len(xs) == 0:
        return False
    if xs[0] != 0:
        return False
    x_max = 0
    for i in range(1, len(xs)):
        if xs[i] > x_max + 1 or xs[i] == xs[i - 1]:
            return False
        elif xs[i] == x_max + 1:
            x_max += 1
    return True


def is_redistribution_of(zs, xs):
    """
    Checks whether zs can be constructed by performing a single redistribution
    of two gametes in xs.
    """
    a, b = None, None
    d = {gz: None for gz in set(zs)}
    for gx, gz in zip(xs, zs):
        if d[gz] is None:
            d[gz] = gx
        elif a is None and gx == d[gz]:
            pass
        elif a is None:
            a, b = gz, gx
        elif gz == a and gx != b:
            return False
        elif d[gz] != a and gx != d[gz]:
            return False
        else:
            pass
    return True


def is_redistribution_of_naive(zs, xs):
    source_gametes_of_gz = {
        gz: {gx for (i, gx) in enumerate(xs) if zs[i] == gz} for gz in set(zs)
    }
    s_0 = None
    for s in source_gametes_of_gz.values():
        if len(s) > 2:
            return False
        if s == 2:
            s_0 = s


def is_redistribution_of_multi(zs, xs):
    """
    Checks whether zs can be constructed by performing (multiple) concurrent
    redistributions. This means that the gametes in xs can be used with replacement to create zs
    """
    d = {gz: (None, None) for gz in set(zs)}
    for gx, gz in zip(xs, zs):
        a, b = d[gz]
        if a is None:
            d[gz] = (gx, None)
        elif b is None and gx == a:
            pass
        elif b is None:
            d[gz] = (a, gx)
        elif gx != b:
            return False
        else:
            pass
    return True


def requires_multipoint(zs, xs, g, i):
    """
    Returns True if the creation of gz[..i] (inclusive) requires multipoint crossover.
    """
    swaps = -1
    prev_value = None
    for gz, gx in zip(zs[: i + 1], xs):
        if gz == g and gx != prev_value:
            prev_value = gx
            swaps += 1
        if swaps > 1:
            return True
    return False


def generate_raw_redistributions_no_filter(xs):
    assert is_distribute(xs)
    n_loci = len(xs)
    n_pop = max(xs) + 1

    def bt(gx, gy, zs, i, new_gamete_values, gz_idx_max):
        if i == n_loci:
            yield zs.copy()
        elif zs[i] is not None:
            yield from bt(gx, gy, zs, i + 1, new_gamete_values, gz_idx_max)
        else:
            for j, gz in enumerate(new_gamete_values[: gz_idx_max + 1]):
                # set value
                zs[i] = gz
                # check if set value would require more than 1 crossover
                if requires_multipoint(zs, xs, gz, i):
                    zs[i] = None
                    continue
                # check if set value would create a gamete that is dominated
                # backtrack
                if j == gz_idx_max:
                    yield from bt(gx, gy, zs, i + 1, new_gamete_values, gz_idx_max + 1)
                else:
                    yield from bt(gx, gy, zs, i + 1, new_gamete_values, gz_idx_max)
                zs[i] = None

    for gx in range(n_pop - 1):
        for gy in range(gx + 1, n_pop):
            # fix non-chosen gametes
            zs = [None if g in (gx, gy) else g for g in xs]
            new_gamete_values = [gx, gy] + list(range(n_pop, n_loci))
            # backtrack on remaining
            yield from bt(gx, gy, zs, 0, new_gamete_values, 0)


def generate_raw_redistributions(xs):
    assert is_distribute(xs)
    n_loci = len(xs)
    out = []
    for zs in sorted(generate_raw_redistributions_no_filter(xs)):
        zs = list(zs)
        if all(gx == gz for gx, gz in zip(xs, zs)):
            continue
        if any(dominates_gametewise(zs, ys) for ys in out):
            continue
        if dominates_gametewise(zs, xs):
            continue
        # remove supersequences of zs
        j = 0
        while j < len(out):
            ys = out[j]
            if dominates_gametewise(ys, zs):
                out.pop(j)
            else:
                j += 1
        out.append(zs)
    return out


def remove_subsequence_if_including(xss, zs):
    # assumes xss is already non-dominating
    if any(is_subsequence_of(ys, zs) for ys in xss):
        return
    # remove supersequences of zs
    j = 0
    while j < len(xss):
        ys = xss[j]
        if is_subsequence_of(zs, ys):
            xss.pop(j)
        else:
            j += 1


def dominates_gametewise(zs, ys):
    """
    Returns True if and only if every gamete in `zs` is dominated
    by a gamete in `ys`.
    """
    d = {gz: None for gz in set(zs)}
    for gz, gy in zip(zs, ys):
        if d[gz] is None:
            d[gz] = gy
        elif d[gz] != gy:
            return False
    return True


def dbg(x):
    __import__("pprint").pprint(x)
    return x


def prune_by_dominance_after_simplification(xss):
    # TODO: if we can prove that there is only two this can be simplified further
    choice = [False] * len(xss)
    # Separate children by size
    sizes = {len(xs) for xs in xss}
    children_by_size = {n: [] for n in sizes}
    for i, xs in enumerate(xss):
        children_by_size[len(xs)].append((i, xs))

    for xss_sub in children_by_size.values():
        xss_sub = prune_by_dominance_func(
            xss_sub,
            func=lambda ixs, iys: dominates_gametewise(iys[1], ixs[1]),
        )
        for i, xs in xss_sub:
            choice[i] = True

    return [xs for xs, b in zip(xss, choice) if b]


def prune_by_dominance_func(xs, func):
    """
    Removes dominated values from `xs` using `func` where `xs` is iterable and
    x dominates y if `func(x, y)`. This function mutates `xs`.
    """
    n = len(xs)
    choice = [True] * n
    for i, x in enumerate(xs):
        if not choice[i]:
            continue
        for j in range(i):
            if choice[j]:
                y = xs[j]
                if func(x, y):
                    choice[j] = False
                if func(y, x):
                    choice[i] = False
                    break
    # Return the chosen elements because that's how we did it in Rust
    return [x for x, b in zip(xs, choice) if b]
    # move elements into their correct positions
    i = j = 0
    while j < n:
        if choice[j]:
            xs[i] = xs[j]
            i += 1
        j += 1
    # and drain excess elements
    for _ in range(i, n):
        xs.pop()


def generate_simplified_redistributions_brute_force(xs):
    d = {}
    for zs in generate_raw_redistributions(xs):
        zs_simplified = _simplify_dist_array(zs, len(zs))
        d[str(zs_simplified)] = zs_simplified
    return prune_by_dominance_after_simplification(list(d.values()))
