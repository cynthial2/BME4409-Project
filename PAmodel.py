# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

import functions.compute_Ra as ra
import functions.compute_ub as basal
import functions.bolus_calc as bolus
import pandas as pd
import numpy as np
from scipy.integrate import odeint


def PA_glc_model(y, t, params, params_ins, basal_values, weight, AC, Ra, u=None, ub=None, Idata=None):
    """ ODE model of glucose-insulin system including meals, insulin injections and PA.

    :param list y: model state
    :param int t: time point
    :param dict params: model parameters
    :param dict params_ins: parameters of insulin kinetics model
    :param list basal_values: basal glucose level, Gb [mg/dl], and basal insulin level, Ib [muU/ml]
    :param float weight: body weight [kg]
    :param float AC: accelerometer count [counts/min]
    :param float Ra: rate of glucose appearance from meals [mg/min]
    :param float u: insulin input [muU/min]. Provide u and ub, or Idata
    :param float ub: basal insulin input [muU/min]. Provide u and ub, or Idata
    :param float Idata: insulin concentration [muU/ml]. Provide u and ub, or Idata
    :return: dydt
    """

    x1, x2, I, X, Q1, Q2, Y, Z, rGU, rGP, tPA, PAint, rdepl, th = y
    Gb, Ib = basal_values

    ''' Parameters '''
    # core model parameters:
    p1 = params['p1']
    p2 = params['p2']
    p3 = params['p3']
    p4 = params['p4']
    p5 = params['p5']
    Vg = params['Vg']

    # insulin parameters
    k1 = params_ins['k1']
    k2 = params_ins['k2']
    k3 = params_ins['k3']
    k4 = params_ins['k4']
    Vi = params_ins['Vi']

    # PA-driven insulin sensitivity parameters:
    tau_AC = params['tau_AC']
    b = params['b']
    tau_Z = params['tau_Z']

    # PA-driven GU & GP parameters
    alpha = params['alpha']
    q1 = params['q1']
    q2 = params['q2']
    q3l = params['q3l']
    q4l = params['q4l']
    q3h = params['q3h']
    q4h = params['q4h']
    q5 = params['q5']

    # Glycogen depletion
    beta = params['beta']
    q6 = params['q6']
    adepl = params['adepl']
    bdepl = params['bdepl']

    # transfer function parameters
    aY = params['aY']
    aAC = params['aAC']
    ah = params['ah']
    n1 = params['n1']
    n2 = params['n2']
    tp = params['tp']

    ''' basal glucose '''
    Qb = Gb * Vg

    ''' basal insulin action '''
    Xb = p3 / p2 * Ib

    ''' Transfer functions '''
    # PA / no PA transfer function
    fY = (Y / aY) ** n1 / (1 + (Y / aY) ** n1)
    fAC = (AC / aAC) ** n2 / (1 + (AC / aAC) ** n2)

    # moderate / high intensity transfer function
    fHI = (AC / ah) ** n2 / (1 + (AC / ah) ** n2)
    fp = (th / tp) ** n2 / (1 + (th / tp) ** n2)

    # depletion time
    if tPA != 0:
        t_depl = - adepl * PAint / tPA + bdepl
    else:
        t_depl = bdepl
    # depletion / no depletion transfer function
    ft = (tPA / t_depl) ** n1 / (1 + (tPA / t_depl) ** n1)

    ''' GP parameters '''
    # GP parameters q3 and q4 based on PA intensity
    q3 = (1 - fp) * q3l + fp * q3h
    q4 = (1 - fp) * q4l + fp * q4h

    # maximum drop in GP after depletion sets in
    rm = beta * (q3 / q4 * Y + (1 - alpha) * (p1 + Xb))

    ''' ODE model '''
    # x1, x2, Ic, X
    if u is not None:
        ins_action = [- k1 * x1 + u + ub,
                      k1 * x1 - (k2 + k3) * x2,
                      k2 / (Vi * weight) * x2 - k4 * I,
                      - p2 * X + p3 * I]
    elif Idata is not None:
        ins_action = [0,
                      0,
                      0,
                      - p2 * X + p3 * Idata]

    # Q1, Q2, Y, Z, rGU, rGP, tPA, PAint, rdepl, th
    glc_dyn = [- (p1 + rGU - (rGP - rdepl) + (1 + Z) * X) * Q1 - p4 * Q1 + p5 * Q2 +
               (p1 + Xb) * Qb + Ra / weight,
               p4 * Q1 - p5 * Q2,

               - 1 / tau_AC * Y + 1 / tau_AC * AC,
               b * fY * Y - (1 - fY) / tau_Z * Z,

               q1 * fY * Y - q2 * rGU,
               q3 * fY * Y - q4 * rGP,

               fAC - (1 - fAC) * tPA,
               fAC * AC - (1 - fAC) * PAint,
               q6 * (ft * rm - rdepl),

               fHI - (1 - fHI) * q5 * th
               ]

    dydt = ins_action + glc_dyn

    return dydt


def compute_glucose_data(params, basal_values, weight, AC, meal, u=None, ub=None, I=None,
                         params_meal={'f': [1], 'tau_m': [40]}):
    """ Compute PA model output for calibration and validation studies.
    Insulin injections (u and ub) or plasma insulin concentration (I) can be given as input.

    :param dict params: model parameters
    :param list basal_values: basal glucose level, Gb [mg/dl], and basal insulin level, Ib [muU/ml]
    :param float weight: body weight [kg]
    :param ndarray AC: vector of accelerometer counts [counts/min]
    :param ndarray meal: vector of meal input [mgCHO/min]
    :param ndarray u: vector of insulin input [muU/min]. Provide u and ub, or I
    :param ndarray ub: vector of basal insulin input [muU/min]. Provide u and ub, or I
    :param ndarray I: vector of insulin concentration [muU/ml]. Provide u and ub, or I
    :param dict params_meal: meal parameters f and tau_m
    :return: model
    """

    # compute basal glucose mass and basal insulin action
    Gb, Ib = basal_values
    Qb = Gb * params['Vg']
    Xb = params['p3'] / params['p2'] * Ib

    # compute meal rate of appearance
    t_meal = np.where(meal != 0)[0]
    D = meal[t_meal]
    Ra = ra.compute_Ra(len(AC), D, t_meal, params_meal['f'], params_meal['tau_m'])

    # load insulin parameters
    params_ins = pd.read_csv('parameters/params_insulin.csv', header=None,
                             dtype={0: str}, delimiter=';').set_index(0).squeeze().to_dict()

    # define y0
    if I is not None:
        y0 = [0, 0, 0, Xb, Qb, params['p4'] / params['p5'] * Qb] + [0] * 8
    elif u is not None:
        y0 = [ub[0]/params_ins['k1'], ub[0]/(params_ins['k2']+params_ins['k3']), Ib, Xb, Qb,
              params['p4'] / params['p5'] * Qb] + [0] * 8

    # initialise model
    model = np.zeros((len(AC), len(y0)))
    model[0, :] = y0

    # if u is given, compute model output including insulin kinetics
    if u is not None:
        for step in range(len(u) - 1):
            t_int = [0, 1]
            sol = odeint(PA_glc_model, y0, t_int,
                         args=(params, params_ins, [Gb, Ib], weight, AC[step], Ra[step], u[step], ub[step], None))
            model[step + 1, :] = sol[1]
            y0 = sol[1]
        model = pd.DataFrame(data=model, columns=['x1', 'x2', 'I', 'X', 'Q1', 'Q2',
                                                  'Y', 'Z', 'rGU', 'rGP', 'tPA', 'PAint', 'rdepl', 'th'])

    # if I is given, compute model output and use insulin profile as input to determine insulin action
    elif I is not None:
        for step in range(len(I) - 1):
            t_int = [0, 1]
            sol = odeint(PA_glc_model, y0, t_int,
                         args=(params, params_ins, [Gb, Ib], weight, AC[step], Ra[step], None, None, I[step]))
            model[step + 1, :] = sol[1]
            y0 = sol[1]
        model = pd.DataFrame(data=model[:, 3:], columns=['X', 'Q1', 'Q2', 'Y', 'Z', 'rGU', 'rGP', 'tPA',
                                                         'PAint', 'rdepl', 'th'])
        model['I'] = I

    # append glucose concentration to dataframe
    model['G'] = model['Q1'] / params['Vg']

    return model


def compute_glucose_pred_interval(params_pred, basal_values, weight, AC, meal, u=None, ub=None, I=None,
                                  params_meal={'f': [1], 'tau_m': [40]}):
    """ Compute prediction intervals for validation studies.
        Insulin injections (u and ub) or plasma insulin concentration (I) can be given as input.

        :param DataFrame params_pred: parameter sets for predictions
        :param list basal_values: basal glucose level, Gb [mg/dl], and basal insulin level, Ib [muU/ml]
        :param float weight: body weight [kg]
        :param ndarray AC: vector of accelerometer counts [counts/min]
        :param ndarray meal: vector of meal input [mgCHO/min]
        :param ndarray u: vector of insulin input [muU/min]. Provide u and ub, or I
        :param ndarray ub: vector of basal insulin input [muU/min]. Provide u and ub, or I
        :param ndarray I: vector of insulin concentration [muU/ml]. Provide u and ub, or I
        :param dict params_meal: meal parameters f and tau_m
        :return: glucose predictions
        """

    G_pred = np.zeros((params_pred.shape[0], len(AC)))

    for i in range(params_pred.shape[0]):
        params_tmp = params_pred.iloc[i].to_dict()

        if len(params_meal['f']) > 1:
            f_tmp = params_meal['f'][i]
            tau_tmp = params_meal['tau_m'][i]
        else:
            f_tmp = params_meal['f']
            tau_tmp = params_meal['tau_m']

        model = compute_glucose_data(params_tmp, basal_values, weight, AC, meal, u=u, ub=ub, I=I,
                                     params_meal={'f': f_tmp, 'tau_m': tau_tmp})
        G_pred[i, :] = np.array(model['G'])

    return G_pred


def compute_glucose_simulation(params, basal_values, weight, AC, meal, u, ub=None,
                               params_meal={'f': [1], 'tau_m': [40]},
                               bolus_calc=False, params_bolus={'Gt': 120, 'ICR': 17, 'CF': 28}):
    """ Compute PA model output for simulation studies.

    :param dict params: model parameters
    :param list basal_values: basal glucose level, Gb [mg/dl], and basal insulin level, Ib [muU/ml]
    :param float weight: body weight [kg]
    :param ndarray AC: vector of accelerometer counts [counts/min]
    :param ndarray meal: vector of meal input [mgCHO/min]
    :param ndarray u: vector of insulin input [muU/min]. If bolus_calc is True, fill with zeros
    :param ndarray ub: vector of basal insulin input [muU/min]. If None, ub is computed from basal insulin level
    :param dict params_meal: meal parameters f and tau_m
    :param bool bolus_calc: If True, meal insulin boluses are computed
    :param dict params_bolus: Parameters for bolus calculation. Target glucose Gt [mg/dl],
                              insulin-to-carb ratio ICR and correction factor CF
    :return: model
    """

    # compute basal glucose mass and basal insulin action
    Gb, Ib = basal_values
    Qb = Gb * params['Vg']
    Xb = params['p3'] / params['p2'] * Ib

    # compute meal rate of appearance
    t_meal = np.where(meal != 0)[0]
    D = meal[t_meal]
    Ra = ra.compute_Ra(len(AC), D, t_meal, params_meal['f'], params_meal['tau_m'])

    # load insulin parameters
    params_ins = pd.read_csv('parameters/params_insulin.csv', header=None,
                             dtype={0: str}, delimiter=';').set_index(0).squeeze().to_dict()

    # compute basal insulin infusion rate if not given
    if ub is None:
        ub0 = basal.compute_ub(Ib, weight, params_ins)
        ub = ub0 * np.ones(len(AC))

    # define y0
    y0 = [ub[0]/params_ins['k1'], ub[0]/(params_ins['k2']+params_ins['k3']), Ib, Xb, Qb,
          params['p4'] / params['p5'] * Qb] + [0] * 8

    # initialise model
    model = np.zeros((len(AC), len(y0)))
    model[0, :] = y0

    # compute model output
    for step in range(len(u) - 1):

        # compute insulin bolus at meal times
        if bolus_calc:
            if step in t_meal:
                u_bol = bolus.bolus_calc(meal[step]*1e-3, model[step, 4]/params['Vg'], params_bolus['Gt'],
                                         params_bolus['ICR'], params_bolus['CF'])
                u[step] += u_bol

        t_int = [0, 1]
        sol = odeint(PA_glc_model, y0, t_int,
                     args=(params, params_ins, [Gb, Ib], weight, AC[step], Ra[step], u[step], ub[step], None))
        model[step + 1, :] = sol[1]
        y0 = sol[1]

    model = pd.DataFrame(data=model, columns=['x1', 'x2', 'I', 'X', 'Q1', 'Q2',
                                              'Y', 'Z', 'rGU', 'rGP', 'tPA', 'PAint', 'rdepl', 'th'])
    model['G'] = model['Q1'] / params['Vg']

    return model


def compute_glucose_replay(params, basal_values, weight, AC, meal, u, ub, G0=None,
                           params_meal={'f': [1], 'tau_m': [40]}):
    """ Compute PA model output for replay simulations.

    :param dict params: patient parameters
    :param list basal_values: basal glucose level, Gb [mg/dl], and basal insulin level, Ib [muU/ml]
    :param float weight: body weight [kg]
    :param ndarray AC: vector of accelerometer counts [counts/min]
    :param ndarray meal: vector of meal input [mgCHO/min]
    :param ndarray u: vector of insulin input [muU/min]
    :param ndarray ub: vector of basal insulin input [muU/min]
    :param float G0: initial glucose level [mg/dl]
    :param dict params_meal: meal parameters f and tau_m
    :return: model
    """

    # compute basal glucose mass and basal insulin action
    Gb, Ib = basal_values
    Qb = Gb * params['Vg']
    Xb = params['p3'] / params['p2'] * Ib

    # compute meal rate of appearance
    t_meal = np.where(meal != 0)[0]
    D = meal[t_meal]
    Ra = ra.compute_Ra(len(AC), D, t_meal, params_meal['f'], params_meal['tau_m'])

    # define initial glucose
    if G0 is None:
        Q0 = Qb
    else:
        Q0 = G0 * params['Vg']

    # load insulin parameters
    params_ins = pd.read_csv('parameters/params_insulin.csv', header=None,
                             dtype={0: str}, delimiter=';').set_index(0).squeeze().to_dict()

    # define y0
    y0 = [ub[0]/params_ins['k1'], ub[0]/(params_ins['k2']+params_ins['k3']), Ib, Xb, Q0,
          params['p4'] / params['p5'] * Q0] + [0] * 8

    # initialise model
    model = np.zeros((len(AC), len(y0)))
    model[0, :] = y0

    for step in range(len(u) - 1):
        t_int = [0, 1]
        sol = odeint(PA_glc_model, y0, t_int,
                     args=(params, params_ins, [Gb, Ib], weight, AC[step], Ra[step], u[step], ub[step], None))
        model[step + 1, :] = sol[1]
        y0 = sol[1]

    model = pd.DataFrame(data=model, columns=['x1', 'x2', 'I', 'X', 'Q1', 'Q2',
                                              'Y', 'Z', 'rGU', 'rGP', 'tPA', 'PAint', 'rdepl', 'th'])
    model['G'] = model['Q1'] / params['Vg']

    return model
