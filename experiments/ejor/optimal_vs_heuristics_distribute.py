"""
In this experiment, we will explore the gap between optimal optimal methods
and the wedges heuristic on randomly generated distribute instances.
"""

import json

import eugene.solvers.base_min_crossings_wedges as ew
import eugene.solvers.base_min_crossings_minizinc as em
import eugene.db as ed

with ed.CachedDB("./distribute_data_2.json", ed.DistributeDB,) as db:

    solver_opt = ed.DBSolver(
        lambda instance: em.breeding_program_distribute(
            len(instance), instance), db
    )

    def solver_heu(instance): return ew.breeding_program(
        len(instance), instance)

    # for n_loci in range(2, 10):
    #     print()
    #     cases = [eu.random_distribute_instance(n_loci) for _ in range(100)]
    #     obj_opt = [solver_opt(arr).objective for arr in cases]
    #     obj_heu = [solver_heu(arr).objective for arr in cases]
    #     for case, opt, heu in zip(cases, obj_opt, obj_heu):
    #         print(case, opt, heu, heu - opt, sep='\t')

    with open("./data/optgap_instances.json") as f:
        case_dict = json.load(f)

    for k, cases in case_dict.items():
        n_loci = str(k)
        obj_arr_opt = []
        obj_arr_heu = []

        for case in cases:
            obj_opt = solver_opt(case).objective
            obj_heu = solver_heu(case).objective
            obj_arr_opt.append(obj_opt)
            obj_arr_heu.append(obj_heu)
            print(case, obj_opt, obj_heu, obj_heu - obj_opt, sep="\t")
