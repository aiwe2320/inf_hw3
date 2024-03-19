from scipy.stats import nbinom

def NB_draws(R0, k, num_draws):
    # Calculate relevant parameters
    mean = R0
    variance = mean + (mean**2)/k
    p = mean/variance
    n = mean**2 / (variance - mean)
    
    # Make draws from negative binomial
    draws = nbinom.rvs(n=n,p=p,size=num_draws)
    
    return draws
    
def run_NB_branching_process(R0, k, I0, Nt_cutoff, num_trials=10):
    '''
    Run branching process simulation
    If number of infections exceeds Nt_cutoff, is considered infinite
    
    Return 1 if infinite, 0 if stops
    '''
    # Store results of simulations
    results = []
    
    # Run a number of trials 
    for x in range(num_trials):
        # Initial draw for I0 individuals
        next_gen = NB_draws(R0, k, I0)
        Nt = I0 + sum(next_gen)
        result = -1
        
        # Loop until process ends or exceeds Nt_cutoff
        while (sum(next_gen) > 0 and Nt < Nt_cutoff):
            nprev = sum(next_gen)  # Store count in previous generation
            next_gen = NB_draws(R0, k, num_draws=nprev)  # Draw next generation
            Nt += sum(next_gen)
        
        # Check outcome
        if (sum(next_gen) == 0):  # Epidemic stops
            result = 0
        elif (Nt > Nt_cutoff):  # Epidemic went forever
            result = 1
        
        # Store result
        results.append(result)
    
    # Return list of results
    return results

def calc_epi_end_prob(results):
    # Count up all times the epidemic stopped
    end_count = 0
    for result in results:
        if (result == 0):
            end_count += 1
    
    q = end_count / len(results)
    return q