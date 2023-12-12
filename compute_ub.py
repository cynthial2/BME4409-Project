# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

import pandas as pd

def compute_ub(Ib, weight, params_ins=None):
    """ Compute basal insulin infusion rate from basal insulin concentration.
    :param float Ib: basal insulin [muU/ml]
    :param float weight: body weight [kg]
    :param dict params_ins: parameters of insulin kinetics model.
    :return: basal insulin infusion rate [muU/min]
    """

    if params_ins is None:
        params_ins = pd.read_csv('parameters/params_insulin.csv', header=None,
                                 dtype={0: str}, delimiter=';').set_index(0).squeeze().to_dict()

    ub = Ib * (params_ins['k2'] + params_ins['k3']) / params_ins['k2'] * (params_ins['k4'] * params_ins['Vi'] * weight)

    return ub
