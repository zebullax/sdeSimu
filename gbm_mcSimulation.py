import argparse
import random
import math 
import numpy
import matplotlib.pylab as plt

# Simulate dynamics for SDE described as 
# dSt = St[(m)dt + (sigma)dWt] , with solution
# St=S0 exp{[m-(sigma^2)/2]t + (sigma)Wt}
def simuGBM(simuParams) :     
    s0      = simuParams["gi"]
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    sigma   = simuParams["gv"]
    mu      = simuParams["gd"] - (sigma**2)/2
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
def simuPoissonProcessJumpTime(simuParams) : 
    t           = 0
    dt          = simuParams["dt"]
    l           = simuParams["pj"]
    T           = simuParams["T"]
    jmpTimes    = []
    F           = 1 - math.exp(-l*dt)
        
    while t <= T :         
        if (random.uniform(0,1) <= F) :
            jmpTimes.append(t)
        t = t + dt
    return jmpTimes

def simuPoissonProcess(simuParams) :    
    results     = [[],[]]
    jumpTimes   = simuPoissonProcessJumpTime(simuParams)

    for jmpTime in jumpTimes :
            results[0].append(jmpTime)
            jumpSize = 1 if simuParams["pu"] else numpy.random.poisson(simuParams["ps"])
            results[1].append(jumpSize)
    return results

def outputSimuToFile(filename, simuResult) :
    if (args["simuGbm"]) :
        outFile = open(filename, mode="w")
        for path in simuResult :
            for val in path : 
                print(val, sep=",", file=outFile)
            print("\n", file=outFile)
    elif (args["simuPoisson"]) :
         for path in simuResult :
            for val in path : 
                print(val, sep=",", file=outFile)

def outputSimuGbmToPlot(simuParams, simuResult) : 
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    nbSteps = int(T/dt + 1) # + 1 for x0
    tScale  = numpy.linspace(0, T, nbSteps)

    mainFrm = plt.figure()
    mainFig = mainFrm.add_subplot(111)
    mainFig.set_title("Geometric Brownian motion process")
    mainFig.set_xlabel("time")
    mainFig.set_ylabel("S(t)")
    for path in simuResult :
        xScale  = numpy.asarray(path)   
        mainFig.plot(tScale, xScale)    
    mainFig.set_xbound(0.0, T*1.2)
    plt.show()
     
def outputSimuPoissonToPlot(simuParams, simuResult) : 
    T       = simuParams["T"]
    dt      = simuParams["dt"]    
    
    mainFrm = plt.figure()
    mainFig = mainFrm.add_subplot(111)
    mainFig.set_title("Poisson process")
    mainFig.set_xlabel("time")
    mainFig.set_ylabel("N(t)")
    mainFig.plot(*simuResult, marker='_', color='k', ls='', markersize=8)    
    mainFig.set_xbound(0.0, T*1.2)    
    jmpMax = max(jmpSz for jmpSz in simuResult[1])+0.2
    mainFig.set_ybound(0.0, jmpMax)
    plt.show()

def main() :
    argparser = argparse.ArgumentParser(prog = "Monte Carlo simulation of GBM",\
       add_help = True)
    #Common params
    argparser.add_argument("--T", "--maturity", type=float, default=1,\
        help="End time for the simulation, unit is left to user consideration."+
            "Eg: T=2, dt=0.5, assuming that unit is [hour] then a timestep is half an hour;"+
            "T=1, dt=1/365, assuming unit is [year] means that a timestep is a day long..."+
            "All coefficients are expressed in the same time unit; eg using the first example"+
            "and sigma=0.2 has the meaning of a volatility = 0.2/hour")
    argparser.add_argument("--N", "--nbPaths", type=int, default=10,\
        help="Number of paths to simulate")
    argparser.add_argument("--dt", "--timestep", type=float, default=1/365,\
        help="Timestep used for time discretization.")
    #GBM Params
    argparser.add_argument("--simuGbm", type=bool, default=True,\
        help="Simulate geometric Brownian motion")    
    argparser.add_argument("--gv", "--vol", type=float, default=0.2,\
        help="Volatility coefficient")
    argparser.add_argument("--gd", "--drift", type=float, default=0.05,\
        help="Drift coefficient)")
    argparser.add_argument("--gi", "--init", type=float, default=100.0,\
        help="Initial value for the process")
    #Poisson Params
    argparser.add_argument("--simuPoisson", type=bool, default=False,\
        help="Simulate Poisson process")
    argparser.add_argument("--pj", "--intensityJumpEvent", type=float, default=10.0,\
        help="Poisson intensity for jump events")
    argparser.add_argument("--pu", "--unitJumpSize", type=bool, default=False,\
        help="True: unit jump size. False: random jump size")
    argparser.add_argument("--ps", "--intensityJumpSize", type=float, default=3.0,\
        help="Average jump size")
    #Output params
    argparser.add_argument("--of", "--outFile", type=str, \
        help="CSV filename for simulation results")
    argparser.add_argument("--og", "--outGraph", type=bool, default=True,\
        help="Plot simulation results")
   
    args = vars(argparser.parse_args())
    if args["simuPoisson"] :  
        simuResult  = simuPoissonProcess(args)        
    elif args["simuGbm"] :
        simuResult  = simuGBM(args)


    if args["of"] :
        outputSimuToFile(args["of"],simuResult)

    if args["og"] :
        if args["simuGbm"] : 
            outputSimuGbmToPlot(args, simuResult)
        elif args["simuPoisson"] :            
            outputSimuPoissonToPlot(args, simuResult)
            
if __name__ == "__main__" :
        main()