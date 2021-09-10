#Aleksandar Doknic July 2021
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
from os import path
from FF_GenOpt_extern import read_gaussian
from FF_GenOpt_extern import read_charmm
from FF_GenOpt_extern import RunMD
from FF_GenOpt_extern import Compute

rng = np.random.default_rng(seed=42)
startTime = time.time()

#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, numOfFocusPoints, samplesPerDimension, maxdims

if(len(sys.argv) > 1):
    configFile = sys.argv[1]
    print("Config file:",configFile)
    mdinp = sys.argv[2]
    mdout = sys.argv[3]
    qmout = sys.argv[4]
    paramfilename = sys.argv[5]
    outputDir = sys.argv[6]
    numOfFocusPoints = int(sys.argv[7])
    samplesPerDimension = int(sys.argv[8])
    maxdims = None
    if(len(sys.argv) > 9):
        maxdims = int(sys.argv[9]) #plot only maxdims dimensions
        
else:
    print("Not enough parameters")

mdexec = "charmm"
np.set_printoptions(precision=12)

CONFIG_CONTENT = ""
PARAMETER_SETTINGS = []

if (path.exists(configFile)):
    print("Config file", configFile ,"found. Reading:")
    with open(configFile,'r') as f:
        lines = f.readlines()
        for line in lines:
            e = " ".join(line.split()).split("#")[0].strip().split(" ")
            if(e[0] == "MDEXEC"):
                mdexec = e[1]
            if(e[0] == "PARAMETER"):
                PARAMETER_SETTINGS.append([e[1],e[2],float(e[3]),float(e[4]),float(e[5])])
            if(e[0] != ''):
                CONFIG_CONTENT += (str(e)+"\n")
                print(e)

else:
    print("Config file does not exist, using default settings.")

class GenOptimizer:
    def __init__(self):
        self.parameterNames = []
        self.parameterType = []
        self.parameterMinVal = np.array([])
        self.parameterMaxVal = np.array([])
        self.parameterInitVal = np.array([])
        self.parameterCurVal = np.array([])
        self.dim = 0
        self.population = []
        #self.failCounter = 0
        self.qmReferenceData = read_gaussian(qmout)
        self.bestFitnessLog = []
        self.bonds = 0
        self.angles = 0
        self.dihedrals = 0
        self.negativeDihedrals = False

    def addParameter(self,name,ptype,pmin,pmax,init_val):
        self.parameterNames.append(name)
        self.parameterType.append(ptype)
        self.parameterMinVal = np.append(self.parameterMinVal,pmin)
        self.parameterMaxVal = np.append(self.parameterMaxVal,pmax)
        self.parameterInitVal = np.append(self.parameterInitVal,init_val)
        self.parameterCurVal = np.append(self.parameterCurVal,init_val)
        self.dim += 1
        if(ptype=="bond"):
            self.bonds += 1
        if(ptype=="angle"):
            self.angles += 1
        if(ptype=="dihedral"):
            self.dihedrals += 1
            if(pmin<0):
                self.negativeDihedrals = True
        
    def normalizeParameterValues(self, plist):
        reslist = []
        for i in range(self.dim):
            reslist.append((plist[i] - self.parameterMinVal[i]) / (self.parameterMaxVal[i] - self.parameterMinVal[i]))
        return np.array(reslist)
    
    def denormalizeParameterValues(self, plist):
        reslist = []
        for i in range(self.dim):
            reslist.append(plist[i] * (self.parameterMaxVal[i] - self.parameterMinVal[i]) + self.parameterMinVal[i])
        return np.array(reslist)
        
    def testFitness(self, plist):
        self.writeParameters(plist)
        try:
            RunMD(mdexec, mdinp, mdout)
            mdfreq, mdX, mdY, mdZ = read_charmm(mdout)
            qmfreq, qmX, qmY, qmZ = self.qmReferenceData
            return Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ)[0]
        except:
            return 2000.0
            #return 9999999.0
            
    def writeParameters(self, plist):
        assert(len(plist)==self.dim)
        output = ""
        for i in range(self.dim):
            output += "set " + self.parameterNames[i] + " " + str(plist[i]) + "\n"
        with open(paramfilename,'w') as f:
            f.write(output)

    def readParameterSettings(self):
        for p in PARAMETER_SETTINGS:
            genopt.addParameter(p[0],p[1],p[2],p[3],p[4])

def createSlices(fn, numOfFocusPoints, samplesPerDimension, dim, sliceRange, maxdims=None):
    if(maxdims == None):
        maxdims = dim
    slices = []
    for s in range(numOfFocusPoints):
        print("Focus point",s+1,"/",numOfFocusPoints)
        fp = np.random.rand(dim)
        for i in range(dim):
            fp[i] = fp[i] * (sliceRange[i][1]-sliceRange[i][0]) #- (sliceRange[i][1]-sliceRange[i][0])/2

        slices.append([])
        for d in range(maxdims):
            print("Dim",d)
            fpc = np.copy(fp)
            slices[s].append([[],[]])
            for i in np.linspace(sliceRange[d][0],sliceRange[d][1],samplesPerDimension):
                #print("i",i)
                fpc[d] = i
                slices[s][d][0].append(i)
                slices[s][d][1].append(fn(fpc))
    return slices

def plotSlices(slices, numOfFocusPoints, dim, maxdims = None):
    if(maxdims == None):
        maxdims = dim
    for d in range(maxdims): #MODIFIED TO 2, should be dim for full slicing
        for s in range(numOfFocusPoints):
            plt.title("dim="+str(d))
            plt.plot(slices[s][d][0],slices[s][d][1],color="grey",alpha=0.5)
        plt.ylim(0,2000)
        plt.savefig(outputDir+"slice_dim"+str(d))
        plt.clf()

def fnsqr(X):
    return np.sum(X**2)

genopt = GenOptimizer()
genopt.readParameterSettings()
sliceRange = []
for i in range(genopt.dim):
    sliceRange.append([genopt.parameterMinVal[i],genopt.parameterMaxVal[i]])
slices = createSlices(genopt.testFitness,numOfFocusPoints=numOfFocusPoints,samplesPerDimension=samplesPerDimension,sliceRange=sliceRange,dim=genopt.dim,maxdims=maxdims)
plotSlices(slices,numOfFocusPoints=numOfFocusPoints,dim=genopt.dim,maxdims=maxdims)

with open(outputDir+"slice_dim.data","w") as out:
    out.write(str(slices))
