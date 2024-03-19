import sir_comp_utils as sircomp
import numpy as np

sv_0 = [0.24975, 0.24975, 0.24975, 0.24975]
iv_0 = [0.00025, 0.00025, 0.00025, 0.00025]
rv_0 = [0, 0, 0, 0]

pv = [1,2,3,4]
w_inv = [4,4,4,4]
gammav = [3,3,3,3]
c_bar = 0.45
C = np.ones((4,4)) * c_bar

dt = 0.1
tmax = 10

t_list, s_lists, i_lists, r_lists = sircomp.sir_compartmental(sv_0, iv_0, rv_0, pv, C, w_inv, gammav, dt, tmax)

xl = 'Time'
yl = 'Proportion of population'
title = 'i(t) Compartments vs Time'
outfile = 'Q2c_i_plot.png'

sircomp.plot_infections_compartmental(t_list=t_list, i_lists=i_lists,
                                      x_label=xl, y_label=yl,
                                      title=title, outfile=outfile)