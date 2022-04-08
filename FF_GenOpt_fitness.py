#Aleksandar Doknic
#2021
import numpy as np
from FF_GenOpt_extern import read_gaussian
from FF_GenOpt_extern import read_charmm
from FF_GenOpt_extern import RunMD
from FF_GenOpt_extern import Compute
from FF_GenOpt_extern import create_context
from FF_GenOpt_extern import update_context
from FF_GenOpt_extern import to_change
from FF_GenOpt_extern import get_varnames
from FF_GenOpt_extern import normal_mode
from pathlib import Path

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
    def __init__(self,paramSettings,extern,psf,crd,params,varfile,mdexec,mdinp,mdout,qmout,paramfilename):
        self.parameterNames = paramSettings.names
        self.dim = paramSettings.dim
        self.numOfEvaluations = 0
        self.numOfGeneticEvaluations = 0
        self.qmReferenceData = read_gaussian(qmout)
        self.psf = Path(psf)
        self.crd = Path(crd)
        self.params = Path(params)
        self.mdexec = mdexec
        self.mdinp = mdinp
        self.mdout = mdout
        self.qmout = qmout
        self.paramfilename = paramfilename
        self.extern = extern
        self.context = None
        self.topology = None
        self.system = None
        self.integrator = None
        self.positions = None
        self.psf = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        if self.extern == "False":
            self.context,self.topology,self.system,self.integrator,self.positions,self.psf = create_context(
                psf,crd,params)
            self.bonds, self.angles, self.dhis = to_change(self.parameterNames, varfile, self.psf)
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
        if self.extern == "True":
            try:
                RunMD(self.mdexec, self.mdinp, self.mdout)
                mdfreq, mdX, mdY, mdZ = read_charmm(self.mdout)
                qmfreq, qmX, qmY, qmZ = self.qmReferenceData
                fitness = Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ)[0]
                return fitness
            except:
                return 9999999.0
        else:
            #try:
            varnames = get_varnames(self.paramfilename)
            update_context(self.system,self.context,varnames,self.bonds,self.angles,self.dihedrals)
            mdfreq, mdX, mdY, mdZ = normal_mode(self.topology, self.system, self.integrator, self.positions)
            #raise StopIteration
            qmfreq, qmX, qmY, qmZ = self.qmReferenceData
            fitness = Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ)[0]
            return fitness
            #except:
            #    return 9999999.0
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
