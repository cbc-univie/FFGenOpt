#Aleksandar Doknic
#2021
import sys
from os import path
import numpy as np

class Config:
    def load(self):
        if(len(sys.argv) > 1):
            configFile = sys.argv[1]
            print("Config file:",configFile)
        else:
            configFile = "FF_GenOpt.cfg"
            print("No config file given. Default file will be used:",configFile)
            configFile = "FF_GenOpt.cfg"
            
        #Default general settings
        self.mdexec = "charmm"
        self.mdinp = "oac.inp" #CHARMM input file that reads from parameters.str
        self.mdout = "opt.log" #CHARMM computer logfile with vibration vectors
        self.qmout = "qmout.log" #QM computed logfile with vibration vector
        self.paramfilename = "parameters.str"
        self.populationfilename = "lastPopulation.txt"
        self.logfilename = "convergence.log"
        np.set_printoptions(precision=12)

        #Some default algorithm settings
        self.GENERATIONS = 100
        self.POPULATION_SIZE = 25
        self.MUTATIONS_PER_GENERATION = 10
        self.STEP_SIZE = 0.05
        self.CROSSOVER_BLX = 0.5
        self.CROSSOVER_SBX = 0.5
        self.CROSSOVER_UNIFORM = 0.0
        self.MINIMUM_IMPROVEMENT = 0.001
        self.MINIMUM_DIVERSITY = 0.0001
        
        self.LOGGING = 0
        self.CREATE_INITIAL_POPULATION = 1
        self.RECOMPUTE_FITNESS = 1
        
        self.CONFIG_CONTENT = ""
        self.PARAMETER_SETTINGS = []

        if (path.exists(configFile)):
            print("Config file", configFile ,"found. Reading:")
            with open(configFile,'r') as f:
                lines = f.readlines()
                for line in lines:
                    e = " ".join(line.split()).split("#")[0].strip().split(" ")
                    if(e[0] == "MDEXEC"):
                        self.mdexec = e[1]
                    if(e[0] == "MDINP"):
                        self.mdinp = e[1]
                    if(e[0] == "MDOUT"):
                        self.mdout = e[1]
                    if(e[0] == "QMOUT"):
                        self.qmout = e[1]
                    if(e[0] == "PARAM_FILENAME"):
                        self.paramfilename = e[1]
                    if(e[0] == "POPULATION_FILENAME"):
                        self.populationfilename = e[1]
                    if(e[0] == "LOG_FILENAME"):
                        self.logfilename = e[1]
                        
                    if(e[0] == "GENERATIONS"):
                        self.GENERATIONS = int(e[1])
                    if(e[0] == "POPULATION_SIZE"):
                        self.POPULATION_SIZE = int(e[1])
                    if(e[0] == "MUTATIONS_PER_GENERATION"):
                        self.MUTATIONS_PER_GENERATION = int(e[1])
                    if(e[0] == "STEP_SIZE"):
                        self.STEP_SIZE = float(e[1])
                    if(e[0] == "CROSSOVER_BLX"):
                        self.CROSSOVER_BLX = float(e[1])
                    if(e[0] == "CROSSOVER_SBX"):
                        self.CROSSOVER_SBX = float(e[1])
                    if(e[0] == "CROSSOVER_UNIFORM"):
                        self.CROSSOVER_UNIFORM = float(e[1])
                    if(e[0] == "MINIMUM_IMPROVEMENT"):
                        self.MINIMUM_IMPROVEMENT = float(e[1])
                    if(e[0] == "MINIMUM_DIVERSITY"):
                        self.MINIMUM_DIVERSITY = float(e[1])
                        
                    if(e[0] == "LOGGING"):
                        self.LOGGING = int(e[1])
                    if(e[0] == "CREATE_INITIAL_POPULATION"):
                        self.CREATE_INITIAL_POPULATION = int(e[1])
                    if(e[0] == "RECOMPUTE_FITNESS"):
                        self.RECOMPUTE_FITNESS = int(e[1])
                    if(e[0] == "PARAMETER"):
                        self.PARAMETER_SETTINGS.append([e[1],e[2],float(e[3]),float(e[4]),float(e[5])])
                    if(e[0] != ''):
                        self.CONFIG_CONTENT += (str(e)+"\n")
                        #print(e)
        return self.CONFIG_CONTENT
