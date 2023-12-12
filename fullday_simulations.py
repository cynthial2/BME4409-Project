# **********************************************************************************************************************
#  Copyright © 2022 ETH Zurich, Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai,
#  Hans-Michael Kaltenbach; D-BSSE; CSB Group
#  All rights reserved. This program and the accompanying materials are made available under the terms of the BSD-3
#  Clause License which accompanies this distribution, and is available at
#  https://gitlab.com/csb.ethz/t1d-exercise-model/-/blob/main/LICENSE
# **********************************************************************************************************************

# This script is used to generate glucose data in full-day simulations using a model of glucose-insulin regulation
# including exercise, meal intake and insulin injections.
# It generates the data for Figure 4 of the following manuscript:
#
# Title:   New model of glucose-insulin regulation characterizes effects of physical activity and facilitates
#          personalized treatment evaluation in children and adults with type 1 diabetes
# Authors: Julia Deichmann, Sara Bachmann, Marie-Anne Burckhardt, Marc Pfister, Gabor Szinnai, Hans-Michael Kaltenbach*
# *Corresponding author:
#          michael.kaltenbach@bsse.ethz.ch
#
# Date:    March 17, 2022
# Author:  Julia Deichmann <julia.deichmann@bsse.ethz.ch>

import functions.PAmodel as GIM
import functions.plot_fct as plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

''' Define basal values '''

Gb = 120                                    # basal glucose level [mg/dl] (normal = 120)
Ib = 12                                     # basal insulin level [muU/ml]
weight = 70                                 # weight [kg]


''' Define scenarios '''

intensity = [0, 2095, 4317, 6539]           # intensity in AC count (rest, 30%, 60%, 90% VO2max)

duration = [300]                            # PA duration [min] (5 hour marathon, ends at 2pm)
pa_start = [240]                            # timing of PA start [min] (start marathon at 9am)

t_meal = [60, 540, 780]                     # timing of meals [min] (baseline, stress, insulin injections)

D = [15e4, 84e3, 84e3]                      # CHO (carbohydrates only) amount [mg] (baseline, stress, insulin injections)

insulin_injection = False

''' Select Stimulation to test'''
print('Menu:\n'
      '1) Baseline\n'
      '2) Insulin Injections at Mealtimes\n'
      '3) Epinephrine Injection (Stress)\n'
      '4) Periodic Carbohydrate Consumption during Marathon\n')

choose = input('Choose a number: ')
choose = int(choose)

if choose == 1:
    pass
elif choose == 2:
    insulin_injection = True
elif choose == 3:
    Gb = 140
elif choose == 4:
    t_meal = [60, 285, 330, 375, 420, 465, 510, 540, 780]  # timing of meals and electrolytes during marathon [min]
    D = [15e4, 22e3, 22e3, 22e3, 22e3, 22e3, 22e3, 84e3,
         84e3]  # CHO (carbohydrates only) amount for meals and electrolytes [mg]
else:
    print('Invalid Choice!')
    exit()


''' Parameters '''

params = pd.read_csv('parameters/params_standard.csv', header=None,
                     dtype={0: str}, delimiter=';').set_index(0).squeeze().to_dict()


''' Define model input '''

t_sim = 1441                               # simulation duration [min]
t = np.arange(t_sim)                        # create time points

meal = np.zeros(t_sim)                      # define meal input [mg/min]
meal[t_meal] = D


''' Compute glucose profile '''

for i in range(len(pa_start)):
    for j in range(len(duration)):
        G = []
        for k in range(len(intensity)):
            if (j == 0) or (k in [0, 1, 2]):

                ''' Insulin injection '''
                if insulin_injection:
                    t_insulin = [60, 420, 780]   # timing of insulin [min] (same as food times)
                    U_amount = [1.28e7, 1.28e7, 1.28e7]
                    u = np.zeros(t_sim)         # define insulin input [muU/min] (amount injected at that minute)
                    u[t_insulin] = U_amount

                ''' No insulin injection '''
                if not insulin_injection:
                    u = np.zeros(t_sim)  # define insulin input [muU/min]

                AC = np.zeros(t_sim)        # define AC input [counts/min]
                AC[pa_start[i]:pa_start[i]+duration[j]] = intensity[k]

                model = GIM.compute_glucose_simulation(params, [Gb, Ib], weight, AC, meal, u=u,
                                                       bolus_calc=True)
                G.append(model['G'])


        # print(G) #final value in each section (at time 1440) is the end glucose value

        plot.plot_simulation(G, labels=['rest', '30%', '60%', '90%'],
                             colors=[plt.cm.PuBu(0.85)] + list(plt.cm.Reds(np.linspace(0.3, 0.9, 3))),
                             pa=[pa_start[i], pa_start[i]+duration[j]])


        ''' Shows mealtimes bars on graph '''
        figure_title = ""
        if choose == 1:
            figure_title = "Baseline"
            plt.vlines([t_meal], ymin=-11.749248318987501, ymax=362.63075179805816, color='r', linestyles='dashed',
                       linewidth=0.5)  # baseline
        elif choose == 2:
            figure_title = "Insulin Injection at Mealtimes"
            plt.vlines([t_meal], ymin=-20.23460935937436, ymax=321.1948332395255, color='r', linestyles='dashed',
                       linewidth=0.5) # insulin injection
        elif choose == 3:
            figure_title = "Epinephrine Injection"
            plt.vlines([t_meal], ymin=-10.064541911488764, ymax=376.03603462276163, color='r', linestyles='dashed',
                       linewidth=0.5)  # stress
        elif choose == 4:
            figure_title = "Consumed Energy Gel during Marathon"
            plt.vlines([t_meal], ymin=-10.064541911488764, ymax=376.03603462276163, color='r', linestyles='dashed',
                       linewidth=0.5)  # energy gel

        ''' Plot visuals'''
        plt.title(figure_title)
        plt.tight_layout()
        plt.savefig('graphs/' + figure_title)
        plt.show()
        print(figure_title, "graph generated")
