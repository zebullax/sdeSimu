import argparse
import random
import math 
import numpy
import matplotlib.pylab as plt

# Simulate dynamics for SDE described as 
# dSt = St[(m)dt + (sigma)dWt] , with solution
# St=S0 exp{[m-(sigma^2)/2]t + (sigma)Wt}
def simuGBM(simuParams) :     
    s0      = simuParams["i"]
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    sigma   = simuParams["s"]
    mu      = simuParams["m"] - (sigma**2)/2
    nbPaths = simuParams["N"]
    nbSteps = int(T/dt)
    results = []
    
    for j in range(nbPaths) :
        simuVal = [s0]
        for i in range(nbSteps) : 
            lastVal = simuVal[len(simuVal)-1]
            simuVal.append(lastVal * math.exp(mu*dt + sigma*math.sqrt(dt)*random.normalvariate(0,sigma)))
        results.append(simuVal[:])
    return results

# Simulate jump times for poisson process
# with intensity l
def simuPoissonProcess(simuParams) : 
    t       = 0
    dt      = simuParams["dt"]
    l       = simuParams["l"]
    T       = simuParams["T"]
    jumps   = []
    F       = 1 - math.exp(-l*dt)
        
    while t <= T :         
        if (random.uniform(0,1) <= F) :
            jumps.append(t)
        t = t + dt
    return jumps

def outputSimuToFile(filename, simuResult) : 
    outFile = open(filename, mode="w")    
    for val in simuResult : 
        print(val, sep=",", file=outFile)

def outputSimuGbmToPlot(simuParams, simuResult) : 
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    nbSteps = int(T/dt + 1) # + 1 for x0
    tScale  = numpy.linspace(0, T, nbSteps)

    for path in simuResult :
        xScale  = numpy.asarray(path)   
        plt.plot(tScale, xScale)
    plt.xlabel("time(y)")
    plt.ylabel("S(t)")
    plt.show()
     
def outputSimuPoissonToPlot(simuParams, simuResult) : 
    T       = simuParams["T"]
    dt      = simuParams["dt"]         
    jumpT   = numpy.asarray(simuResult)
    
    for time in jumpT : 
        plt.plot(time,1,"o")
    plt.xlabel("time(y)")
    plt.ylabel("N(t)")
    plt.show()

def main() :
    argparser = argparse.ArgumentParser(prog = "Monte Carlo simulation of GBM", add_help = True)
    #Common params
    argparser.add_argument("--T", "--maturity", type=float, default=1.0,\
        help="End time (unit year) for the simulation")
    argparser.add_argument("--N", "--nbPaths", type=int, default=100,\
        help="Number of paths to simulate")
    #GBM Params
    argparser.add_argument("--gbm", type=bool, default=True,\
        help="Simulate geometric Brownian motion")    
    argparser.add_argument("--s", "--std", type=float, default=0.2,\
        help="Annual standard deviation (coefficient for the dWt term)")
    argparser.add_argument("--m", "--mean", type=float, default=0.0,\
        help="Annual mean (coefficient for the dt term)")
    argparser.add_argument("--i", "--init", type=float, default=100.0,\
        help="Initial value for the process")    
    argparser.add_argument("--dt", "--timestep", type=float, default=0.01,\
        help="Timestep used for time discretization")
    #Poisson Params
    argparser.add_argument("--poisson", type=bool, default=False,\
        help="Simulate Poisson process")
    argparser.add_argument("--l", "--lambda", type=float, default=10.0,\
        help="Jump intensity for Poisson process")
    #Output params
    argparser.add_argument("--o", "--outFile", type=str, \
        help="CSV filename for simulation results")
    argparser.add_argument("--g", "--outGraph", type=bool, default=True,\
        help="Plot simulation results")
   
    args = vars(argparser.parse_args())
    if args["poisson"] : 
        simuResult  = simuPoissonProcess(args)
    elif args["gbm"] :
        simuResult  = simuGBM(args)


    if args["o"] :
        outputSimuToFile(args["o"],simuResult)
    if args["g"] :
        if args["gbm"] : 
            outputSimuGbmToPlot(args, simuResult)
        elif args["poisson"] : 
            outputSimuPoissonToPlot(args, simuResult)
            
if __name__ == "__main__" :
        main()