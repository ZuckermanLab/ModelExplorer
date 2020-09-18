# Latin Hypercube Sampling - August George - 4/17/2020

from smt.sampling_methods import LHS
import numpy as np


def latin_hypercube(x_limits, n):
    #  see smt library tutorial online:
    #  https://smt.readthedocs.io/en/latest/_src_docs/sampling_methods/lhs.html
    sampling = LHS(xlimits=x_limits, criterion='c')  # create hypercube and select from center of cube partitions (using package)
    x = sampling(n)  #  sample n points in Latin hypercube (using package)
    return x  


def main():
    
    ### example usage
    n = 10  # draw 10 samples from the latin hypercube
    x_lim = np.array([[-10,10],[0,2],[-1,1]])  # parameter ranges for 3 dim
    y = latin_hypercube(x_lim,n)  # generate sample points using Latin hypercube
    print("3D Latin hypercube example:")
    print(y.shape)  # 10 points in a 3D sample space
    print(y)  # samples
 


if __name__ == "__main__": #best practice to use main
    main()