# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

import numpy as np

def bolus_calc(D, BG, BGt, ICR, CF):
    """ Compute meal insulin bolus.
    :param float D: CHO amount [g]
    :param float BG: current glucose level [mg/dl]
    :param float BGt: glucose target [mg/dl]
    :param float ICR: insulin-to-carb ratio
    :param float CF: correction factor
    :return: insulin bolus [muU]
    """

    bolus = D / ICR + (BG - BGt) / CF
    bolus = np.round(2 * bolus) / 2 * 1e6
    bolus = max(bolus, 0)

    return bolus
