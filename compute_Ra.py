# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

import numpy as np

def compute_Ra(dur, D, t_meal, f, tau_m):
    """ Compute glucose rate of appearance from meal CHOs.
    :param int dur: simulation duration [min]
    :param list D: CHO amount of meals [mg]
    :param list t_meal: time of meal intake [min]
    :param list f: bioavailability. One shared value for all meals, or one value per meal.
    :param list tau_m: time of maximum glucose appearance [min]. One shared value for all meals, or one value per meal.
    :return: rate of glucose appearance [mg/min]
    """

    n_meal = len(D)
    Ra = np.zeros(dur)

    for i in range(n_meal):
        dur_meal = dur - t_meal[i]
        t = np.arange(dur_meal)

        if len(f) == 1:
            Ra[t_meal[i]:] += f[0] * D[i] * t * np.exp(-t / tau_m[0]) / tau_m[0] ** 2
        else:
            Ra[t_meal[i]:] += f[i] * D[i] * t * np.exp(-t / tau_m[i]) / tau_m[i] ** 2

    return Ra
