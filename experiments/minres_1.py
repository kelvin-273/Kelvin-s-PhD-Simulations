"""
In this experiment we are going to investigate the solution quality and runtime
required by the following algorithms:
- Segment algorithm
- Distribute A*
- MinRes MIP
- MinRes Sampling
- MinRes Greedy

We will test these algorithms on random and distribute instances (Distribute A*
will be omitted from testing on random instances).

How will we test them
"""

import eugene_pywrapper
import eugene.utils as eu
import eugene.solvers.base_min_resources_mip as emip
import eugene.plant_models.plant2 as ep2
import pymongo as pm
import argparse
from random import random, seed
from collections import namedtuple
import time
from multiprocessing import Process
from multiprocessing import Pipe

TIMEOUT = 300
THREAD_DELTA = 0.125

LOCI = list(range(2, 10))
ALLOWED_LOCI = tuple(LOCI)
INSTANCE_TYPES = ("random", "distribute")


def run_with_timeout(f, args=(), timeout=None):
    tx, rx = Pipe()
    res = None
    p = Process(group=None, target=f, args=(args, tx))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
    else:
        res = rx.recv()
    p.join()
    return res


def solver_segment_dist_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_gen.segment.breeding_program_distribute_minres_python(
        *args
    )
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_segment_dist(xs, r, gamma):
    return run_with_timeout(
        solver_segment_dist_aux, (xs, None), TIMEOUT + THREAD_DELTA
    )


def solver_distastar_dist_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_cross.distribute_astar.breeding_program_distribute_minres_python(
        *args
    )
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_distastar_dist(xs, r, gamma):
    return run_with_timeout(
        solver_distastar_dist_aux, (xs, None), TIMEOUT + THREAD_DELTA
    )


def solver_greedy_dist_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_res.greedy_dom.breeding_program_distribute_python(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_greedy_dist(xs, r, gamma):
    return run_with_timeout(
        solver_greedy_dist_aux, (xs, r, gamma, None), TIMEOUT + THREAD_DELTA
    )


def solver_sampling_dist_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_res.sampling.breeding_program_distribute_python(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_sampling_dist(xs, r, gamma):
    return run_with_timeout(
        solver_sampling_dist_aux, (xs, r, gamma, None), TIMEOUT + THREAD_DELTA
    )


def solver_mip_dist_aux(args, tx):
    start = time.time()
    cs = emip.breeding_program(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_mip_dist(xs, r, gamma):
    n_loci = xs[0].n_loci
    return run_with_timeout(
        solver_mip_dist_aux, (n_loci, xs, r, gamma), TIMEOUT + THREAD_DELTA
    )


def solver_segment_rand_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_gen.segment.breeding_program_minres_python(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_segment_rand(n_loci, pop_0, r, gamma):
    return run_with_timeout(
        solver_segment_rand_aux, (n_loci, pop_0, None), TIMEOUT + THREAD_DELTA
    )


def solver_greedy_rand_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_res.greedy_dom.breeding_program_python(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_greedy_rand(n_loci, pop_0, r, gamma):
    return run_with_timeout(
        solver_greedy_rand_aux, (n_loci, pop_0, r, gamma, None), TIMEOUT + THREAD_DELTA
    )


def solver_sampling_rand_aux(args, tx):
    start = time.time()
    cs = eugene_pywrapper.min_res.sampling.breeding_program_python(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_sampling_rand(n_loci, pop_0, r, gamma):
    return run_with_timeout(
        solver_sampling_rand_aux,
        (n_loci, pop_0, r, gamma, None),
        TIMEOUT + THREAD_DELTA,
    )


def solver_mip_rand_aux(args, tx):
    start = time.time()
    cs = emip.breeding_program(*args)
    res_time = time.time() - start
    tx.send((cs, res_time))
    tx.close()


def solver_mip_rand(n_loci, pop_0, r, gamma):
    return run_with_timeout(
        solver_mip_rand_aux, (n_loci, pop_0, r, gamma), TIMEOUT + THREAD_DELTA
    )


def build_instances():
    seed(0)
    instances_random = {
        n_loci: [
            (
                ep2.PlantSPCBitarray.initial_pop_random(
                    n_loci=n_loci, n_individuals=n_individuals
                ),
                [random() / 2 for _ in range(n_loci - 1)],
                0.5 + random() / 2,
            )
            for _ in range(100)
            for n_individuals in [2, 4, 6, 8]
        ]
        for n_loci in LOCI
    }
    instances_distribute = {
        n_loci: [
            (
                eu.random_distribute_instance(n_loci),
                [random() / 2 for _ in range(n_loci - 1)],
                0.5 + random() / 2,
            )
            for _ in range(100)
        ]
        for n_loci in LOCI
    }
    return instances_distribute, instances_random


def get_results_collection(reset_collection=False):
    client = pm.MongoClient("mongodb://localhost:27017/")
    db = client["eugene_breeding_programs"]
    if reset_collection:
        db.drop_collection("minres_experiments_1")
    return db["minres_experiments_1"]


# print the contents of the collection
def print_collection_contents(collection):
    for document in collection.find():
        print(document)


Solver = namedtuple(
    "Solver",
    [
        "name",
        "f_rand",
        "f_dist",
        "f_obj",
        "f_rand_pop_conversion",
        "f_dist_pop_conversion",
    ],
)
SOLVERS = [
    Solver(
        "Distribute A*",
        None,
        solver_distastar_dist,
        lambda solution, r, gamma: solution.resources(r, gamma),
        None,
        None,
    ),
    Solver(
        "Segment Algorithm",
        solver_segment_rand,
        solver_segment_dist,
        lambda solution, r, gamma: solution.resources(r, gamma),
        lambda pop_0: [x.to_bool_list() for x in pop_0],
        None,
    ),
    Solver(
        "MinRes Greedy",
        solver_greedy_rand,
        solver_greedy_dist,
        lambda solution, r, gamma: solution.resources(r, gamma),
        lambda pop_0: [x.to_bool_list() for x in pop_0],
        None,
    ),
    Solver(
        "MinRes Sampling",
        solver_sampling_rand,
        solver_sampling_dist,
        lambda solution, r, gamma: solution.resources(r, gamma),
        lambda pop_0: [x.to_bool_list() for x in pop_0],
        None,
    ),
    Solver(
        "MinRes MIP",
        solver_mip_rand,
        solver_mip_dist,
        lambda solution, r, gamma: solution.objective,
        lambda pop_0: [x.to_PlantSPC() for x in pop_0],
        eu.distribute_to_plants,
    ),
]


def select_solvers(solver_names):
    if not solver_names:
        return SOLVERS

    requested_names = []
    for solver_name in solver_names:
        requested_names.extend(
            name.strip() for name in solver_name.split(",") if name.strip()
        )

    solvers_by_name = {solver.name: solver for solver in SOLVERS}
    unknown_names = sorted(set(requested_names) - set(solvers_by_name))
    if unknown_names:
        available_names = ", ".join(solvers_by_name)
        unknown_names = ", ".join(unknown_names)
        raise ValueError(
            f"Unknown solver(s): {unknown_names}. Available solvers: {available_names}"
        )

    return [solvers_by_name[name] for name in requested_names]


def select_instance_types(instance_types):
    if not instance_types:
        return ["random"]

    requested_types = []
    for instance_type in instance_types:
        requested_types.extend(
            name.strip() for name in instance_type.split(",") if name.strip()
        )

    unknown_types = sorted(set(requested_types) - set(INSTANCE_TYPES))
    if unknown_types:
        available_types = ", ".join(INSTANCE_TYPES)
        unknown_types = ", ".join(unknown_types)
        raise ValueError(
            f"Unknown instance type(s): {unknown_types}. "
            f"Available instance types: {available_types}"
        )

    return requested_types


def select_loci(loci_args):
    if not loci_args:
        return LOCI

    selected_loci = []
    for loci_arg in loci_args:
        for token in loci_arg.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                start_text, end_text = token.split("-", 1)
                start = int(start_text)
                end = int(end_text)
                if start > end:
                    raise ValueError(
                        f"Invalid loci range: {token}. Range start must be <= end."
                    )
                selected_loci.extend(range(start, end + 1))
            else:
                selected_loci.append(int(token))

    invalid_loci = sorted(set(selected_loci) - set(ALLOWED_LOCI))
    if invalid_loci:
        allowed_text = f"{ALLOWED_LOCI[0]}-{ALLOWED_LOCI[-1]}"
        invalid_text = ", ".join(str(value) for value in invalid_loci)
        raise ValueError(
            f"Unsupported loci value(s): {invalid_text}. "
            f"Allowed loci are within {allowed_text}."
        )

    return sorted(set(selected_loci))


def result_exists(collection, solver, n_loci, instance_type, instance_number):
    return collection.find_one(
        {
            "solver": solver.name,
            "n_loci": n_loci,
            "instance_type": instance_type,
            "instance_number": instance_number,
        }
    )


def store_result(
    collection, solver, n_loci, instance_type, instance_number, objective, res_time
):
    collection.insert_one(
        {
            "solver": solver.name,
            "n_loci": n_loci,
            "instance_type": instance_type,
            "instance_number": instance_number,
            "objective": objective,
            "time": res_time,
        }
    )


def run_random_instances(collection, solvers, instances_random, selected_loci):
    for n_loci in selected_loci:
        for solver in solvers:
            for i, (pop_0, rec_rates, gamma) in enumerate(instances_random[n_loci]):
                print(
                    f"Running {solver.name} on random instance with n_loci={n_loci}, instance={i}"
                )
                if solver.f_rand is None:
                    print(
                        f"Skipping {solver.name} on random instance with n_loci={n_loci}, instance={i} (no random solver)"
                    )
                    continue
                if result_exists(collection, solver, n_loci, "random", i):
                    print(
                        f"Skipping {solver.name} on random instance with n_loci={n_loci}, instance={i} (already exists)"
                    )
                    continue
                pop = solver.f_rand_pop_conversion(pop_0)
                result = solver.f_rand(n_loci, pop, rec_rates, gamma)
                if result is None:
                    objective = None
                    res_time = TIMEOUT
                else:
                    solution, res_time = result
                    objective = solver.f_obj(solution, rec_rates, gamma)
                print(
                    f"Result for {solver.name} on random instance with n_loci={n_loci}, instance={i}: objective={objective}, time={res_time}"
                )
                store_result(
                    collection, solver, n_loci, "random", i, objective, res_time
                )


def run_distribute_instances(
    collection, solvers, instances_distribute, selected_loci
):
    for n_loci in selected_loci:
        for solver in solvers:
            for i, (dist_array, rec_rates, gamma) in enumerate(
                instances_distribute[n_loci]
            ):
                print(
                    f"Running {solver.name} on distribute instance with n_loci={n_loci}, instance={i}"
                )
                if solver.f_dist is None:
                    print(
                        f"Skipping {solver.name} on distribute instance with n_loci={n_loci}, instance={i} (no distribute solver)"
                    )
                    continue
                if result_exists(collection, solver, n_loci, "distribute", i):
                    print(
                        f"Skipping {solver.name} on distribute instance with n_loci={n_loci}, instance={i} (already exists)"
                    )
                    continue
                instance = (
                    solver.f_dist_pop_conversion(dist_array)
                    if solver.f_dist_pop_conversion is not None
                    else dist_array
                )
                result = solver.f_dist(instance, rec_rates, gamma)
                if result is None:
                    objective = None
                    res_time = TIMEOUT
                else:
                    solution, res_time = result
                    objective = solver.f_obj(solution, rec_rates, gamma)
                print(
                    f"Result for {solver.name} on distribute instance with n_loci={n_loci}, instance={i}: objective={objective}, time={res_time}"
                )
                store_result(
                    collection, solver, n_loci, "distribute", i, objective, res_time
                )


def main(reset_collection=False, solver_names=None, instance_types=None, loci_args=None):
    solvers = select_solvers(solver_names)
    selected_instance_types = select_instance_types(instance_types)
    selected_loci = select_loci(loci_args)
    instances_distribute, instances_random = build_instances()
    collection = get_results_collection(reset_collection)
    print_collection_contents(collection)

    if "random" in selected_instance_types:
        run_random_instances(collection, solvers, instances_random, selected_loci)
    if "distribute" in selected_instance_types:
        run_distribute_instances(
            collection, solvers, instances_distribute, selected_loci
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reset-collection",
        action="store_true",
        help="Drop minres_experiments_1 before starting the benchmark.",
    )
    parser.add_argument(
        "--solvers",
        nargs="+",
        metavar="SOLVER",
        help=(
            "Only run the selected solver names. Accepts quoted names and "
            "comma-separated lists. Available: "
            + ", ".join(solver.name for solver in SOLVERS)
        ),
    )
    parser.add_argument(
        "--instance-types",
        nargs="+",
        metavar="TYPE",
        help=(
            "Only run the selected instance types. Accepts random, distribute, "
            "quoted names, and comma-separated lists. Defaults to random."
        ),
    )
    parser.add_argument(
        "--loci",
        nargs="+",
        metavar="LOCUS",
        help=(
            "Only run the selected loci counts. Accepts values and inclusive "
            f"ranges like 2 4 6, 2-5, or 2,4,6-8. Allowed loci: {ALLOWED_LOCI[0]}-{ALLOWED_LOCI[-1]}."
        ),
    )
    args = parser.parse_args()
    try:
        main(
            reset_collection=args.reset_collection,
            solver_names=args.solvers,
            instance_types=args.instance_types,
            loci_args=args.loci,
        )
    except ValueError as exc:
        parser.error(str(exc))
