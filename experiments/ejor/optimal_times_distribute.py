"""
In this experiment, we ask how different solvers perform on randomly generated
distribute instances.

We'll measure the expected the run-time for randomly generated distribute
instances. As the number of distribute instances for n_loci = 1..6 is less
than 100, we shall simply take the expected value over all loci. The number
instances where n_loci = 7..8 is less than 1000, so these can be chosen from
the full set. Instances where n_loci ≥ 9can sampled by generating and checking
against a list of previously observed instances.
"""

import os
import json
from random import seed
from multiprocessing import Process, Pipe, set_start_method

import eugene.utils as eu
import eugene.solvers.base_min_crossings_astar as east
import eugene.solvers.base_min_crossings_mip as emip
import eugene.solvers.base_min_crossings_minizinc as emzn
import eugene.solvers.base_min_crossings_distribute_astar as edis

seed(0)

N_LOCI = list(range(2, 11))
N_INST = 100
TIMEOUT = 300
THREAD_DELTA = 0.125

INSTANCES = {
    n_loci: [eu.random_distribute_instance(n_loci) for _ in range(N_INST)]
    for n_loci in N_LOCI
}


# INSTANCES = {int(k): v for k, v in json.load(open("./data/optgap_instances.json")).items()}


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


def solver_distastar_aux(args, tx):
    start = time.time()
    res_obj = edis.breeding_program_distribute(*args)
    res_time = time.time() - start
    tx.send((res_obj, res_time))
    tx.close()


def solver_distastar(n_loci, pop_0):
    return run_with_timeout(
        solver_distastar_aux, (pop_0,), TIMEOUT + THREAD_DELTA
    )


def solver_astar_aux(args, tx):
    start = time.time()
    res_obj = east.breeding_program_distribute(*args)
    res_time = time.time() - start
    tx.send((res_obj, res_time))
    tx.close()


def solver_astar(n_loci, pop_0):
    return run_with_timeout(
        solver_astar_aux, (pop_0,), TIMEOUT + THREAD_DELTA
    )


CTX_SAT = emzn.MinizincContext.from_solver_and_model_file(
    "cp-sat", "./eugene/solvers/minizinc/mincross.mzn"
)


def solver_cp_aux(args, tx):
    start = time.time()
    res_obj = emzn.breeding_program_distribute(*args)
    res_time = time.time() - start
    tx.send((res_obj, res_time))
    tx.close()


def solver_cp_sat(n_loci, pop_0):
    return run_with_timeout(
        solver_cp_aux, (n_loci, pop_0, CTX_SAT), TIMEOUT + THREAD_DELTA
    )


CTX_MIP = (DeprecationWarning("gurobi is expired"), None)[1]
# emzn.MinizincContext.from_solver_and_model_file(
#     "gurobi", "./eugene/solvers/minizinc/mincross.mzn"
# )


def solver_cp_mip(n_loci, pop_0):
    return run_with_timeout(
        solver_cp_aux, (n_loci, pop_0, CTX_MIP), TIMEOUT + THREAD_DELTA
    )


def solver_mip_aux(args, tx):
    start = time.time()
    res_obj = emip.breeding_program_distribute(*args)
    res_time = time.time() - start
    tx.send((res_obj, res_time))
    tx.close()


def solver_mip(n_loci, pop_0):
    return run_with_timeout(
        solver_mip_aux, (pop_0,), TIMEOUT + THREAD_DELTA
    )


SOLVERS = {
    # "CANZAR": None,
    "CP-SAT": solver_cp_sat,
    "CP-MIP": solver_cp_mip,
    "MIP": solver_mip,
    "A*": solver_astar,
    "Dist-A*": solver_distastar,
}

if __name__ == "__main__":
    set_start_method("fork")
    import argparse
    import time

    parser = argparse.ArgumentParser(
        "Experiment runner for MinCross algorithms on random instances"
    )
    parser.add_argument("-s", "--solver")
    parser.add_argument("-o", "--output_file")
    parser.add_argument("-l", "--loci")
    args = parser.parse_args()
    print(args)

    solver = args.solver
    output_file = args.output_file
    loci_string = args.loci

    if output_file is None:
        output_file = ""

    if os.path.exists(output_file):
        print(f"output_file exists: {output_file}")
        exit()

    output_dir = os.path.dirname(output_file)
    if output_dir != "" and not os.path.exists(output_dir):
        # make directory
        os.makedirs(output_dir)
        # then we can assume that output_file does not exist

    if solver is None:
        solvers = SOLVERS.items()
    elif solver.upper() == "DIST-ASTAR":
        solvers = [("Dist-A*", SOLVERS["Dist-A*"])]
    elif solver.upper() == "ASTAR":
        solvers = [("A*", SOLVERS["A*"])]
    elif solver.upper() not in SOLVERS:
        raise ValueError("solver not available")
    else:
        solvers = [(solver.upper(), SOLVERS[solver.upper()])]

    if loci_string is None:
        loci_string = "2-8"

    if loci_string.isnumeric():
        N_LOCI = [eval(loci_string)]
    elif loci_string.count("-") == 1 and all(
        s.isnumeric() for s in loci_string.split("-")
    ):
        s, e = (int(c) for c in loci_string.split("-"))
        N_LOCI = list(range(s, e + 1))
    elif loci_string.count("-") == 0 and all(
        s.isnumeric() for s in loci_string.split(",")
    ):
        N_LOCI = [eval(s) for s in loci_string.split(",")]
    else:
        raise ValueError("incorrect format for loci set")

    for name, solver in solvers:
        for n_loci, instances in ((n_loci, INSTANCES[n_loci]) for n_loci in N_LOCI):
            objectives = [None] * len(instances)
            times_l1 = [0] * len(instances)
            start = time.time()
            for i, inst in enumerate(instances):
                result = solver(n_loci, inst)
                if result is None:
                    objectives[i] = None
                    times_l1[i] = TIMEOUT
                else:
                    res_obj, res_time = result
                    if res_time >= TIMEOUT:
                        objectives[i] = None
                        times_l1[i] = TIMEOUT
                    else:
                        objectives[i] = res_obj.objective
                        times_l1[i] = res_time
            if output_file == "":
                print(
                    json.dumps(
                        {
                            "solver": name,
                            "n_loci": n_loci,
                            "times_l1": times_l1,
                            "objectives": objectives,
                        }
                    )
                )
            else:
                with open(output_file, "a") as f:
                    print(
                        json.dumps(
                            {
                                "solver": name,
                                "n_loci": n_loci,
                                "times_l1": times_l1,
                                "objectives": objectives,
                            }
                        ),
                        file=f,
                    )
