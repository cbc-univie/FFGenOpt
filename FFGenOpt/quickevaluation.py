#Aleksandar Doknic
#Computes the fitness score for the input vector
#2021
import numpy as np
import sys
from FF_GenOpt_extern import read_gaussian
from FF_GenOpt_extern import read_charmm
from FF_GenOpt_extern import RunMD
from FF_GenOpt_extern import Compute
from FF_GenOpt_conf import Config

class TestFitnessFunction:
    def __init__(self,paramSettings):
        self.dim = paramSettings.dim
        self.target = np.random.rand(self.dim)*1000
        self.numOfEvaluations = 0
        self.numOfGeneticEvaluations = 0
        #print("Target",self.target)
    def compute(self,parameters,typeOfEvaluation="unknown"):
        assert(type(parameters)==np.ndarray)
        self.numOfEvaluations += 1
        if(typeOfEvaluation=="genetic"):
            self.numOfGeneticEvaluations += 1
        return np.sum((parameters-self.target)**2) #test

class ParameterSettings:
    def __init__(self):
        self.names = []
        self.types = []
        self.minVals = np.array([])
        self.maxVals = np.array([])
        self.initVals = np.array([])
        self.dim = 0
    def addParameterSetting(self,name,ptype,pmin,pmax,initVal):
        self.names.append(name)
        self.types.append(ptype)
        self.minVals = np.append(self.minVals,pmin)
        self.maxVals = np.append(self.maxVals,pmax)
        self.initVals = np.append(self.initVals,initVal)
        self.dim += 1
    def normalizeParameterValues(self, plist):
        reslist = []
        for i in range(self.dim):
            reslist.append((plist[i] - self.minVals[i]) / (self.maxVals[i] - self.minVals[i]))
        return np.array(reslist)
    def denormalizeParameterValues(self, plist):
        reslist = []
        for i in range(self.dim):
            reslist.append(plist[i] * (self.maxVals[i] - self.minVals[i]) + self.minVals[i])
        return np.array(reslist)

class FitnessFunction:
    def __init__(self,paramSettings,mdexec,mdinp,mdout,qmout,paramfilename):
        self.parameterNames = paramSettings.names
        self.dim = paramSettings.dim
        self.numOfEvaluations = 0
        self.numOfGeneticEvaluations = 0
        self.qmReferenceData = read_gaussian(qmout)
        self.mdexec = mdexec
        self.mdinp = mdinp
        self.mdout = mdout
        self.qmout = qmout
        self.paramfilename = paramfilename
        self.paramSettings = paramSettings
    def writeParameters(self, plist):
        assert(len(plist)==self.dim)
        output = ""
        for i in range(self.dim):
            output += "set " + self.parameterNames[i] + " " + str(plist[i]) + "\n"
        with open(self.paramfilename,'w') as f:
            f.write(output)
    def compute(self,parameters,typeOfEvaluation="unknown"):
        self.numOfEvaluations += 1
        if(typeOfEvaluation=="genetic"):
            self.numOfGeneticEvaluations += 1
        self.writeParameters(parameters)
        try:
            RunMD(self.mdexec, self.mdinp, self.mdout)
            mdfreq, mdX, mdY, mdZ = read_charmm(self.mdout)
            qmfreq, qmX, qmY, qmZ = self.qmReferenceData
            fitness = Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ)[0]
            return fitness
        except:
            return 9999999.0
    def dumpParametersCHARMM(self,p):
        output = ""
        for i in range(self.paramSettings.dim):
            output += "set " + self.paramSettings.names[i] + " " + "{:.2f}".format(params[i]) + "\n"
        return output
    def dumpFrequencies(self, parameters):
        output = ""
        self.compute(parameters)
        qfreq = np.array(read_gaussian(self.qmout)[0])*0.957 #qm factor from external file
        mdfreq = np.array(read_charmm(self.mdout)[0])
        #print freq
        for i in range(len(qfreq)):
           output += (str("{:10.4f}".format(qfreq[i])) + "\t:" + str("{:10.4f}".format(mdfreq[i])))
           output += "\n"
        output += "\n"
        #print normalized freq
        #for i in range(len(qfreq)):
        #   output += (str("{:10.4f}".format(mdfreq[i]/qfreq[i]))) #+ "\t:" + str("{:10.4f}".format(mdfreq[i])))
        #   output += "\n"
        return output
        
conf = Config()
conf.load()
paramSettings = ParameterSettings()
for p in conf.PARAMETER_SETTINGS:
    paramSettings.addParameterSetting(name=p[0],ptype=p[1],pmin=p[2],pmax=p[3],initVal=p[4])
fitnessFunction = FitnessFunction(paramSettings,conf.mdexec,conf.mdinp,conf.mdout,conf.qmout,conf.paramfilename)
params = [float(i) for i in sys.argv[2].replace("]","").replace("[","").split(",")]
print("Fitness",fitnessFunction.compute(params),"\n")
print(fitnessFunction.dumpFrequencies(params))
print(fitnessFunction.dumpParametersCHARMM(params))
