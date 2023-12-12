# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
plt.rcParams.update({'font.size': 8})
plt.rcParams['errorbar.capsize'] = 3


def plot_simulation(G, labels, colors=None, linestyles=None, pa=None):
    """ Plot glucose trajectories of 24h simulations.
        G: List of glucose trajectories. """

    t = np.arange(len(G[0]))
    nG = len(G)

    if colors is None:
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, nG))
    if linestyles is None:
        linestyles = nG * [(0, ())]

    xticks_labels = ['06:00', '12:00', '18:00', '00:00']
    x_ticks_pos = np.arange(0, 1440, 6 * 60)

    fig, ax = plt.subplots(1, 1, figsize=(7, 3.4))

    for i in range(nG):
        ax.plot(t, G[i], label=labels[i], linewidth=1.5, color=colors[i], zorder=nG-i, linestyle=linestyles[i])

    ymin, ymax = ax.get_ylim()
    ax.fill_between([t[0], t[-1]], 0, 70, alpha=.2, color='grey')
    ax.fill_between([t[0], t[-1]], 180, 400, alpha=.2, color='grey')

    ''' Shows start/stop of exercise'''
    if pa is not None:
        ax.vlines(pa, ymin-10, ymax+10, color='k', linewidth=1)
        # print(ymin-10, ymax+10)   # check length of red line that indicates CHO intake

    plt.xticks(x_ticks_pos, xticks_labels)
    plt.ylabel('glucose [mg/dl]')
    plt.xlabel('time [hours]')

    plt.xlim(t[0], t[-1])
    plt.ylim(ymin-10, ymax+10)

    plt.legend(framealpha=0)


def plot_personalization(data, model, input):
    """ Plot model output and patient data including CHO, insulin and PA inputs. """

    G_model = model['G']
    I_model = model['I']
    Y_model = model['Y']

    t = np.array(input['time'])
    AC = np.array(input['AC'])
    u = np.array(input['u'])
    ub = np.array(input['ub'])
    meal = np.array(input['meal'])

    idx_c = np.where(meal != 0)[0]
    idx_i = np.where(u != 0)[0]

    xtick_pos = [0, 360, 720, 1080]
    xtick_lab = ['00:00', '06:00', '12:00', '18:00']
    colors = [plt.cm.PuBu(0.85), plt.cm.inferno(0.6)]
    colors2 = plt.cm.viridis([0.1, 0.62])

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(5, 2.8))

    ax[0].plot(t, G_model, zorder=11, c=colors[0], label='model')
    ax[0].plot(data['time'], data['Glc'], zorder=10, c=colors[1], linewidth=1.5, label='data')
    ymin, ymax = (20, 390)
    ax[0].vlines(np.where(input['PA']==1)[0], ymin, ymax, color='grey', alpha=.1)
    ax[0].fill_between([0, 1440], 20, 70, color='grey', alpha=0.2)
    ax[0].fill_between([0, 1440], 180, 390, color='grey', alpha=0.2)
    ax[0].set_ylim(ymin, ymax)
    ax[0].legend(ncol=2, loc='upper left', framealpha=0)
    ax[0].set_ylabel('glucose\n[mg/dl]')

    ax1 = ax[0].twinx()
    ax1.scatter(idx_c, meal[idx_c]*1e-3, label='CHO [g]', color=colors2[1], s=16, marker='x')
    ax1.set_ylabel('CHO [g]')
    ax1.set_ylim(0, 110)
    ax1.set_yticks([0, 50])

    ax[1].plot(t, I_model, zorder=11, c=colors[0], label='model')
    ax[1].set_ylabel('insulin\n[$\mu$U/ml]')
    ax[1].set_ylim(10, 140)

    ax2 = ax[1].twinx()
    ax2.scatter(idx_i, u[idx_i]*1e-6, color=colors2[1], marker='x', s=16, zorder=2)
    ax2.plot(t, ub*1e-6*60, color='k', linewidth=.75, linestyle='--', zorder=1, label='u$_b$ [U/h]')
    ax2.legend(loc='upper left', framealpha=0)
    ax2.set_ylim(0, 11)
    ax2.set_yticks([0, 5, 10])
    ax2.set_ylabel('bolus [U]')

    ax[2].scatter(t, AC, s=0.75, zorder=10, color='k', label='AC')
    ax[2].plot(t, Y_model, c=colors[0], label='Y')
    ymin, ymax = (0, 6000)
    ax[2].vlines(np.where(input['PA']==1)[0], ymin, ymax, color='grey', alpha=.1)
    ax[2].hlines(1500, 0, 1440, color='black', linewidth=.75)
    ax[2].set_xticks(xtick_pos)
    ax[2].set_xticklabels(xtick_lab)
    ax[2].set_ylabel('PA [c/min]')
    ax[2].set_xlim(0, t[-1]+1)
    ax[2].set_ylim(ymin, ymax)
    ax[2].legend(loc='upper left', framealpha=0)


def plot_replay(G, labels, colors=None, linestyles=None, pa=None):
    """ Plot replay simulations.
        G: List of glucose trajectories. """

    t = np.arange(len(G[0]))
    nG = len(G)

    if colors is None:
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, nG))
    if linestyles is None:
        linestyles = nG * [(0, ())]

    xticks_labels = ['00:00', '06:00', '12:00', '18:00']
    x_ticks_pos = np.arange(0, 1440, 6 * 60)

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 1.7))

    for i in range(nG):
        ax.plot(t, G[i], label=labels[i], linewidth=1.5, color=colors[i], linestyle=linestyles[i])

    ymin, ymax = ax.get_ylim()
    ax.fill_between([t[0], t[-1]], 0, 70, alpha=.2, color='grey')
    ax.fill_between([t[0], t[-1]], 180, 400, alpha=.2, color='grey')

    if pa is not None:
        ax.vlines(pa, ymin-10, ymax+10, color='k', linewidth=1)

    plt.xticks(x_ticks_pos, xticks_labels)
    plt.ylabel('glucose [mg/dl]')

    plt.xlim(t[0], t[-1])
    plt.ylim(ymin-10, ymax+10)

    if nG > 4:
        plt.legend(framealpha=0, ncol=2)
    else:
        plt.legend(framealpha=0)


def plot_model(t, model, basal_values, params, pa=None, depl=None, glc_data=None, glc_sd=None, GU=None, GP=None):
    """ Plot PA model output and data including glucose levels, GU and GP. """

    p1 = params['p1']
    p2 = params['p2']
    p3 = params['p3']
    Vg = params['Vg']
    alpha = params['alpha']

    # determine basal levels of insulin action and glucose mass
    glc = model['G']
    Gb, Ib = basal_values
    Qb = Gb * Vg
    Xb = p3 / p2 * Ib

    # compute GU & GP at rest, with insulin-dependent PA effects, and full PA impact
    GU_rest = (1 - alpha) * (p1 + model['X']) * model['Q1']
    GU_PA = (1 - alpha) * (p1 + (1 + model['Z']) * model['X']) * model['Q1'] + model['rGU'] * model['Q1']
    GU_id = GU_rest + (1 - alpha) * model['Z'] * model['X'] * model['Q1']

    GP_rest = (p1 + Xb) * Qb - alpha * (p1 + model['X']) * model['Q1']
    GP_PA = (p1 + Xb) * Qb - alpha * (p1 + (1 + model['Z']) * model['X']) * model['Q1'] + \
            (model['rGP'] - model['rdepl']) * model['Q1']
    GP_id = GP_rest - alpha * model['Z'] * model['X'] * model['Q1']

    # create figure
    colors = [plt.cm.PuBu(0.85), plt.cm.inferno(0.6)]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(3.5, 3.1))

    ax1.plot(t, glc, label='model', linewidth=1.5, c=colors[0])
    if (glc_data is not None) & (glc_sd is not None):
        ax1.errorbar(glc_data['time'], glc_data['Glc'], yerr=glc_sd, c=colors[1],
                 zorder=10, linestyle='', marker='.', markersize=6, linewidth=1, label='data',
                 markeredgewidth=1)
    elif (glc_data is not None) & (glc_sd is None):
        ax1.scatter(glc_data['time'], glc_data['Glc'], c=colors[1], zorder=10, s=12, label='data')
    ymin, ymax = ax1.get_ylim()
    ax1.vlines(pa, ymin-10, ymax+10, color='k', linewidth=.75)
    ax1.set_xlim(t[0]-2, t[-1]+2)
    ax1.set_ylim(ymin-10, ymax+10)
    if depl is not None:
        ax1.vlines(depl, ymin-10, ymax+10, color='k', linestyle=(0, (5, 5)), linewidth=.75)
    ax1.set_ylabel('glucose \n[mg/dl]')
    ax1.legend(framealpha=0, loc=1, ncol=2, bbox_to_anchor=(0.92, 1.05))

    ax2.plot(t, GU_PA, linewidth=1.5, label='PA', zorder=10, c=colors[0])
    if GU is not None:
        ax2.errorbar(glc_data['time'], GU['GU'], yerr=GU['GU_sd'], c=colors[1],
                 zorder=10, linestyle='', marker='.', markersize=6, linewidth=1,
                 markeredgewidth=1)
    ax2.plot(t, GU_rest, linewidth=1.5, label='rest', c=colors[0], linestyle='--')
    ax2.plot(t, GU_id, linewidth=1, color='black')
    ax2.fill_between(t, GU_rest, GU_id, alpha=0, hatch='//', zorder=9)
    ax2.fill_between(t, GU_id, GU_PA, alpha=0.3, color='grey', zorder=8)
    ymin, ymax = ax2.get_ylim()
    ax2.vlines(pa, ymin, ymax, color='black', linewidth=.75)
    ax2.set_ylim(ymin, ymax)
    ax2.set_ylabel('GU \n[mg/kg/min]')
    ax2.legend(framealpha=0)

    ax3.plot(t, GP_PA, linewidth=1.5, zorder=10, c=colors[0])
    if GP is not None:
        ax3.errorbar(glc_data['time'], GP['GP'], yerr=GP['GP_sd'], c=colors[1],
                     zorder=10, linestyle='', marker='.', markersize=6, linewidth=1,
                     markeredgewidth=1)
    ax3.plot(t, GP_rest, linewidth=1.5, color=colors[0], linestyle='--')
    ax3.plot(t, GP_id, linewidth=1, color='black')
    ax3.fill_between(t, GP_rest, GP_id, alpha=0, hatch='//', label='id', zorder=9)
    ax3.fill_between(t, GP_id, GP_PA, alpha=0.3, color='grey', label='ii', zorder=8)
    ymin, ymax = ax3.get_ylim()
    ax3.vlines(pa, ymin, ymax, color='black', linewidth=.75)
    ax3.set_ylim(ymin, ymax)

    ax3.set_ylabel('GP \n[mg/kg/min]')
    ax3.set_xlabel('time [min]')
    ax3.legend(framealpha=0)


def plot_glc(t, G, Gb, pa=None, depl=None, glc_data=None, glc_sd=None, Gpred=None, dG=False):
    """ Plot PA model output and data as glucose concentration, G, or difference to basal, G-Gb.
        Prediction intervals can be included. """

    if dG:
        G = G - Gb
        glc_data = glc_data.rename(columns={'dG': 'Glc'})
        if Gpred is not None:
            Gpred = Gpred - Gb

    colors = [plt.cm.PuBu(0.85), plt.cm.inferno(0.6)]

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 1.7))


    ax.plot(t, G, label='model', linewidth=2, color=colors[0], zorder=3)
    if (glc_data is not None) & (glc_sd is not None):
        ax.errorbar(glc_data['time'], glc_data['Glc'], yerr=glc_sd, color=colors[1],
                    zorder=10, linestyle='', marker='.', markersize=6, linewidth=1, label='data',
                    markeredgewidth=1)
    elif (glc_data is not None) & (glc_sd is None):
        ax.scatter(glc_data['time'], glc_data['Glc'], color=colors[1], zorder=10, s=12, label='data')
    if Gpred is not None:
        ax.fill_between(t, np.percentile(Gpred, 2.5, axis=0), np.percentile(Gpred, 97.5, axis=0),
                        color=colors[0], alpha=.3, zorder=2)

    ymin, ymax = ax.get_ylim()
    if pa is not None:
        ax.vlines(pa, ymin-10, ymax+10, color='k', linewidth=1)
    if depl is not None:
        ax.vlines(depl, ymin-10, ymax+10, color='k', linestyle=(0, (5, 5)), linewidth=.75)
    if dG:
        ax.hlines(0, t[0]-2, t[-1]+2, linestyles='--', color='black')
    ax.set_xlim(t[0]-2, t[-1]+2)
    ax.set_ylim(ymin-10, ymax+10)

    ax.set_xlabel('time [min]')
    if dG:
        ax.set_ylabel('$\Delta$G [mg/dl]')
    else:
        ax.set_ylabel('glucose [mg/dl]')


