import eugene_pywrapper
import pandas as pd
from random import random, randrange
from math import log2

GAMMA = 0.9

if __name__ == "__main__":
    samples = []
    for n_loci in range(2, 10):
        for _ in range(100):
            rates = [random() * 0.5 / (n_loci - 1) for _ in range(n_loci - 1)]
            rr_multi = eugene_pywrapper.utils.PyRecRate(rates)
            rr_single = rr_multi.to_singlepoint()
            # rr_single = eugene_pywrapper.utils.PySinglePointRecProb(rates)
            for _ in range(100):
                x_u = randrange(1 << n_loci)
                x_l = randrange(1 << n_loci)

                y_u = randrange(1 << n_loci)
                y_l = randrange(1 << n_loci)

                x = (x_u << n_loci) | x_l
                y = (y_u << n_loci) | y_l
                # random feasible crossing
                if random() < 0.5:
                    x_u, x_l = x_l, x_u
                if random() < 0.5:
                    y_u, y_l = y_l, y_u

                k_x = randrange(n_loci)
                k_y = randrange(n_loci)

                # leaving it this way as opposed to (1 << (n_loci - k_x))
                # because the sampling is symmetrical and is easier to read
                mask_x_u = (1 << k_x) - 1
                mask_x_l = ((1 << n_loci) - 1) ^ mask_x_u
                mask_y_u = (1 << k_y) - 1
                mask_y_l = ((1 << n_loci) - 1) ^ mask_y_u

                z_u = (mask_x_u & x_u) | (mask_x_l & x_l)
                z_l = (mask_y_u & y_u) | (mask_y_l & y_l)

                z = (z_u << n_loci) | z_l

                resources_multi = rr_multi.cost_of_crossing(GAMMA, x, y, z)
                resources_single = rr_single.cost_of_crossing(GAMMA, x, y, z)
                samples.append(
                    {
                        "n_loci": n_loci,
                        "rates": rates,
                        "x": x,
                        "y": y,
                        "z": z,
                        "resources_multi": resources_multi,
                        "resources_single": resources_single,
                        "ratio": log2(resources_single / resources_multi),
                    }
                )

    df = pd.DataFrame(samples)
    print(df)
