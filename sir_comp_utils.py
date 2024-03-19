import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def generate_t_list(t0, dt, tmax):
    t_list = []
    # If 0<dt<1, can't use range function
    if (dt > 0 and dt < 1):
        count = t0
        while (count < tmax):
            t_list.append(count)
            count += dt
    else:  # Get list of times using range
        t_list = list(range(t0, tmax, dt))
    return t_list

def sir_compartmental(sv_0, iv_0, rv_0, pv, C, w_inv, gammav, dt, tmax):
    '''
    Run a compartmental SIR simulation for a system with population structure
    
    sv_0 : Vector of initial s values
    iv_0 : Vector of initial i values
    rv_0 : Vector of initial r values
    pv  : Vector of transmission probabilities
    C : Contact matrix among groups
    w_inv : Vector of 1 / omegas
    gammav : Vector of gammas
    '''
    
    # Establish dynamic matrix + vector
    D_s = np.diag(sv_0)
    i = np.array(iv_0)
    r = np.array(rv_0)
    
    # Set up initial values
    s_lists = []
    i_lists = []
    r_lists = []
    for ndx, i_val in enumerate(iv_0):
        s_lists.append([sv_0[ndx]])
        i_lists.append([i_val])
        r_lists.append([rv_0[ndx]])
    
    # Set up diagonal matrices of constants
    D_p = np.diag(pv)
    D_w_inv = np.diag(w_inv)
    D_gamma = np.diag(gammav)
    
    # Generate list of times for calculations
    t_list = generate_t_list(t0=0, dt=dt, tmax=tmax)
    
    # Step forward in time by calculating ds and di
    for time in t_list:
        # Derivatives
        A = D_s @ D_p @ C @ D_w_inv
        ds = -(A) @ i
        di = (A) @ i - D_gamma @ i
        dr = D_gamma @ i
    
        # Update values
        for ndx, ds_val in enumerate(ds):
            s_prev = D_s[ndx][ndx]
            s_new = ds_val * dt + s_prev 
            D_s[ndx][ndx] = s_new
            s_lists[ndx].append(s_new)
            
            i_prev = i[ndx]
            i_new = di[ndx] * dt + i_prev 
            i[ndx] = i_new
            i_lists[ndx].append(i_new)
            
            r_prev = r[ndx]
            r_new = dr[ndx] * dt + r_prev 
            r[ndx] = r_new
            r_lists[ndx].append(r_new)
    
    # Complete list of times, return results
    t_list.append(t_list[-1] + dt)  # Accounts for situations where tmax is indivisible by dt
    
    return t_list, s_lists, i_lists, r_lists

def plot_infections_compartmental(t_list, i_lists, x_label, y_label, title, outfile):
    '''
    Plot the different infectious compartments versus time with varying opacity
    '''
    num_groups = len(i_lists)
    # Generate labels for each group
    labels = ['Group ' + str(ndx + 1) for ndx in range(len(i_lists))]
    # Generate an even distribution of opacities for each i line
    # The last i group will be the most solid (alpha=1)
    opacities = [(ndx + 1) * (0.5 / num_groups) for ndx in range(len(i_lists))]
    
    fig, ax = plt.subplots()
    for ndx, i_list in enumerate(i_lists): # Plot i data using a loop
        ax.plot(t_list, i_list, label=labels[ndx], alpha=opacities[ndx], color='b')
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()

    plt.savefig(outfile, bbox_inches='tight')