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

    simuVal = [s0]
    nbSteps = int(T/dt)
    for i in range(nbSteps) : 
        lastVal = simuVal[len(simuVal)-1]
        simuVal.append(lastVal * math.exp(mu*dt + sigma*math.sqrt(dt)*random.normalvariate(0,sigma)))
    return simuVal

def outputSimuToFile(filename, simuResult) : 
    outFile = open(filename, mode="w")    
    for val in simuResult : 
        print(val, sep=",", file=outFile)

def outputSimuToPlot(simuParams, simuResult) : 
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    nbSteps = int(T/dt + 1) # + 1 for x0
    tScale  = numpy.linspace(0, T, nbSteps)
    xScale  = numpy.asarray(simuResult)
    plt.plot(tScale, xScale)
    plt.show()
     
def main() :
    argparser = argparse.ArgumentParser(prog = "Monte Carlo simulation of GBM", add_help = True)
    argparser.add_argument("--s", "--std", type=float, default=0.2,\
        help="Standard deviation coefficient for the dWt term")
    argparser.add_argument("--m", "--mean", type=float, default=0.0,\
        help="Mean coefficient for the dt term")
    argparser.add_argument("--i", "--init", type=float, default=100.0,\
        help="Initial value for the SP")
    argparser.add_argument("--T", "--maturity", type=float, default=1.0,\
        help="Upper time value for the simulation")
    argparser.add_argument("--dt", "--timestep", type=float, default=0.01,\
        help="Timestep value used for time discretization")
    argparser.add_argument("--o", "--outFile", type=str, \
        help="Comma delimited file for simulation output results")
    argparser.add_argument("--g", "--outGraph", type=bool, default=True,\
        help="Comma delimited file for simulation output results")
   
    args        = vars(argparser.parse_args())    
    simuResult  = simuGBM(args)

    if args["o"] :
        outputSimuToFile(args["o"],simuResult)
    if args["g"] :
        outputSimuToPlot(args, simuResult)
            
if __name__ == "__main__" :
        main()