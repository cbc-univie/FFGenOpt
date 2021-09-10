#Aleksandar Doknic
#2021
import numpy as np
from FF_GenOpt_extern import read_gaussian
from FF_GenOpt_extern import read_charmm
from FF_GenOpt_extern import RunMD
from FF_GenOpt_extern import Compute

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
        for i in range(len(qfreq)):
           output += (str("{:10.4f}".format(mdfreq[i]/qfreq[i]))) #+ "\t:" + str("{:10.4f}".format(mdfreq[i])))
           output += "\n"
        return output
