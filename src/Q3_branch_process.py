import branch_process_utils as bp

R0 = 3
k_list = [0.1, 0.5, 1, 5, 10]
Nt_cutoff = 100000
I0 = 1
num_trials = 5000

q_list = []

# Run case study
for k in k_list:
    results = bp.run_NB_branching_process(R0, k, I0, Nt_cutoff, num_trials=num_trials)
    q_temp = bp.calc_epi_end_prob(results)
    
    q_list.append(q_temp)

print(q_list)