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

# Generate SIR simulation
t_list, s_lists, i_lists, r_lists = sircomp.sir_compartmental(sv_0, iv_0, rv_0, pv, C, w_inv, gammav, dt, tmax)

# Calculate average probability of transmission over time
pbar_list = sircomp.calculate_average_probability(pv=pv, t_list=t_list, s_lists=s_lists)

# Make susceptibles plot
xl_s = 'Time'
yl_s = 'Proportion of population'
title_s = 's(t) Compartments vs Time'
outfile_s = 'Q2d_s_plot.png'

sircomp.plot_susceptible_compartmental(t_list=t_list, s_lists=s_lists,
                                       x_label=xl_s, y_label=yl_s,
                                       title=title_s, outfile=outfile_s)

# Make probability plot
xl_p = 'Time'
yl_p = 'Proportion of population'
title_p = 'pbar(t) vs Time'
outfile_p = 'Q2d_p_plot.png'

sircomp.plot_probability(t_list=t_list, pbar_list=pbar_list,
                         x_label=xl_p, y_label=yl_p,
                         title=title_p, outfile=outfile_p)