import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import time
from math import e
start_time = time.time()

def test_acf(runs = 1, scale = 10, size = 10):
    '''
    Runs multiple experiments and calculates statistics (avg and var) for the
    correlation time and exponent term. 
    Parameters:
        runs: int (default = 1)
            How many experiments to run. Default is 1.
        scale: int (default = 10)
            Factor for exponential distribution: 1/scale * exp(-x/scale). 
        size: int (default = 100) 
            How many random numbers to draw from exp. distribution
    Returns:
        out: ndarray 
            Contains runs, scale, size, correlation time avg, 
            correlation time variance, fitted exponential term avg, and
            fitted exponential term variance. 

    Notes: Expect ACF to decay exponentially given the trajectories are from a two 
    state Poisson process (exponential distribution).
    Correlation time is the index (time) of the the ACF before the first negative ACF.
    Fitted exponential term is found by the slope of log(ACF) using polyfit
    '''
    exp_fit = np.zeros(runs)
    exp_fit_r = np.zeros(runs)
    c_time = np.zeros(runs).astype(int)
    for i in range(runs):
        data = generate_data(scale, size)
        ACF = acf(data)
        c_time[i] = (np.where(ACF < 0)[0][0] - 1)  # num before first negative number
        # truncate after correlation time. np.s_ allows for [:] notation
        tr_ACF = np.delete(ACF, np.s_[int(c_time[i]+1):],)
        log_ACF = np.log(tr_ACF)  #
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
            np.arange(log_ACF.size), log_ACF)
        exp_fit[i] = slope
        exp_fit_r[i] = r_value*r_value
    #graph(np.arange(log_ACF.size), log_ACF)
    graph(np.arange(tr_ACF.size),tr_ACF, scale)
    return np.array([runs, scale, size, np.mean(c_time), np.var(c_time), np.mean(exp_fit), np.var(exp_fit), np.mean(exp_fit_r)])


def generate_data(scale = 10, size = 10):
    '''
    Generates synthetic trajectory data from a two state Poisson process. 
    Parameters:
        scale: int (default = 10)
            Factor for exponential distribution: 1/scale * exp(-x/scale).
        size: int (default = 100)
            How many random numbers to draw from exp. distribution
    Returns:
        out: ndarray
            Output data contains a series of 0's and 1's representing the duration at 
            each state. 
    Notes: Time between Poisson process can be modeled using the exponential distribution: 
    1/B exp(-x/B) where x is a random variable and B is the scaling factor. Numbers are 
    randomly picked from the exponential distribution, corresponding to the wait times
    for a state transition, and then appended to the trajectory. For example, 7 is picked 
    for state 0, so 7 0's are appended to the list. Then 3 is picked, so 3 1's are appended. 
    This repeats n = size times. 
    '''
    state = 0
    traj = []
    nums = np.ceil(np.random.exponential(scale, size)).astype(int)
    for i in nums:
        if state == 0:
            for j in range(i):
                traj.append(0)
                state = 1
        else:
            for k in range(i):
                traj.append(1)
                state = 0
    return np.asarray(traj)


def acf(data):
    '''
    Calculates the normalized Autocorrelation function (ACF)
    Parameters:
        data: ndarray
            n x m matrix (n of vectors) of data
    Returns:
        out: ndarray
            ACF of data
    Notes: CF_k = <(x(t)-x.mean).(x(t+k)-x.mean)> / <(x(t)-x.mean).(x(t)-x.mean)>
    where k is the lag time, and < > denotes the average over t. 
    '''
    #N = np.shape(data)[0]
    N = data.size
    mean = np.mean(data)
    acf = np.zeros(N)

    # calculate den = <(x(i)-x.mean).(x(i)-x.mean)>
    den = 0.0
    for t in range(N):  # t = 0 to N-1
        den = den + np.dot((data[t]-mean), (data[t]-mean))
    den = den/N  

    # calculate num = <(x(t)-x.mean).(x(t+k)-x.mean)>
    for k in range(N):  # k = 0 to N-1
        num = 0.0
        for t in range(N-k):
            num = num + np.dot((data[t]-mean), (data[t+k]-mean))
        num = num/(N-k)    
        acf[k] = num/den
    return acf


def graph(x,y,scale=1):
    plt.figure(figsize=(10,5))
    f_x = []

    factor = 1.0/scale
    for i in x:
        f_x.append(1*np.exp(-i*2*factor))
    f_x = np.asarray(f_x)
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(x, y, 'b.', label="synthetic data")
    plt.plot(x, f_x, 'r--', label="exp(-2/%s)" % scale)  # extra 2 factor from sum of both states
    plt.title(
        "ACF(x) for two state Poisson process\nExponential distribution: 1/%s *f(x) = exp(-x/%s)" % (scale, scale))
    plt.ylabel("ACF(x)")
    plt.xlabel("Lag time x (truncated)")
    plt.legend()
    plt.subplot(212)
    plt.plot(x, np.log(y), 'b.', label="log[synthetic data]")
    plt.plot(x, np.log(f_x), 'r--', label="log[exp(-2/%s)]" % scale)
    plt.ylabel("log[ACF(x)]")
    plt.xlabel("Lag time x (truncated)")
    plt.legend()
    plt.tight_layout()
    plt.show()


#t = [10, 100, 1000]

t = [1000]
for i in t:
    results = test_acf(runs = 1, scale=i, size=10)
    print("\n")
    print("Randomly selected %s numbers from the exponential distribution: (1/%s)*exp(-x/%s). Repeated %s times." % (int(results[2]), int(results[1]), int(results[1]), int(results[0])))
    print("Correlation Time: Avg = %s, Var = %s." % (results[3], results[4]))
    print("Fitted Exponential Term: Avg = %s, Var = %s, Avg R^2 = %s." % (results[5], results[6], results[7]))
    print("%s (s)" % (time.time() - start_time))
    print("\n")

