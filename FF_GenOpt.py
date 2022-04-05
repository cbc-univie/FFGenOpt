import matplotlib.pyplot as plt
import numpy as np
import time
import datetime
from os import path
from FF_GenOpt_conf import Config
from FF_GenOpt_fitness import TestFitnessFunction
from FF_GenOpt_fitness import FitnessFunction
np.set_printoptions(precision=2)

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

class PopulationMember:
    def __init__(self,values,fitness):
        assert(type(values) == np.ndarray)
        self.values = values
        self.fitness = fitness

class Population:
    def __init__(self,settings,fn,maxSize,mutationsPerGeneration):
        assert(type(settings) == ParameterSettings)
        self.members = []
        self.settings = settings
        self.fn = fn #fitness function
        self.maxSize = maxSize
        self.mutationsPerGeneration = mutationsPerGeneration
    def addPopulationMember(self,member):
        assert(type(member) == PopulationMember)
        self.members.append(member)
    def sortPopulation(self):
        self.members = sorted(self.members,key=lambda x: x.fitness)
    def dumpPopulation(self):
        out = ""
        for i in range(len(self.members)):
            out += "{:.2f}".format(self.members[i].fitness)
            out += ','
            out += str([round(i,3) for i in self.members[i].values])
            out += '\n'
        return out
    def getFitnessList(self):
        return np.array([m.fitness for m in self.members])
    def getFitnessSum(self):
        return np.sum([m.fitness for m in self.members])
    def getFitnessMean(self):
        return self.getFitnessSum()/len(self.members)
    def getFitnessWorst(self):
        return self.members[:-1].fitness
    def getFitnessBest(self):
        return self.members[0].fitness
    def generateInitialPopulation(self,size=10,type="uniform"):
        for i in range(size):
            randVec = np.array([])
            if(type=="uniform"):
                randVec = np.random.rand(self.settings.dim)
            if(type=="gaussian"):
                randVec = np.abs(np.random.normal(0.5,0.2,self.settings.dim))
            randVec = self.settings.denormalizeParameterValues(randVec)
            fitness = self.fn.compute(randVec,"genetic")
            self.addPopulationMember(PopulationMember(randVec,fitness))
            print("Generating initial population",i+1,"/",size," fitness:",fitness)
    def cutPopulation(self,size=None):
        if(size == None):
            size = self.maxSize
        self.sortPopulation()
        self.members = self.members[:size]
    #Take 3 random dimensions and compare for 20 percent of population
    def estimateDiversity(self,sampleDim=3,sampleSize=0.2):
        numOfSamples = max(1,int(len(self.members)*sampleSize))
        msd = 0.0
        for i in range(sampleDim):
            randDim = np.random.randint(0,self.settings.dim)
            for s in range(numOfSamples):
                randSample1 = np.random.randint(0,numOfSamples)
                randSample2 = np.random.randint(0,numOfSamples)
                while(randSample2==randSample1): #in case it is the same
                    randSample2 = np.random.randint(0,numOfSamples)
                msd += (self.members[randSample1].values[randDim]-self.members[randSample2].values[randDim])**2
        msd /= (numOfSamples*sampleDim)
        return msd
    def readPopulation(self,fileName):
        if(path.exists(fileName)):
            print("Found population file!")
            print(fileName)
            print("Recomputing fitness... this may take a while.")
            with open(fileName,'r') as f:
                lines = f.readlines()[1:]
                for line in lines:
                    if (line.strip()==""): #only take before the freq part
                        break
                    #fitness = float(line.replace("\n","").split(",[")[0])
                    plist = [float(p) for p in line.replace("\n","").split(",[")[1].replace("]","").split(",")]
                    print(plist)
                    fitness = self.fn.compute(plist)
                    self.addPopulationMember(PopulationMember(np.array(plist),fitness))
                    print("Found population member with fitness",fitness)
            print("Total population:",len(self.members))
            return True
        else:
            print("No population file found! Creating a new population.")
            return False
        
class GeneticOperators:
    def __init__(self,paramSettings):
        self.paramSettings = paramSettings
    def uniformSelection(self,population):
        return np.random.choice(population.members)
    def wheelSelection(self,population):
        population.sortPopulation()
        probabilities = population.getFitnessList()
        probabilities = np.flip(probabilities)
        total = np.sum(probabilities)
        probabilities /= total
        return np.random.choice(population.members,p=probabilities)
    def rankSelection(self,population):
        total = len(population.members)
        probabilities = [i for i in range(total,0,-1)]
        probabilities /= np.sum(probabilities)
        return np.random.choice(population.members,p=probabilities)
    def crossoverUniform(self,p1,p2):
        assert(len(p1)==len(p2))
        o = []
        for i in range(len(p1)):
            if(i%2==0):
                o.append(p1[i])
            else:
                o.append(p2[i])
        return np.array(o)
    def crossoverBLX(self,p1,p2,alpha=0.5):
        assert(len(p1)==len(p2))
        o = []
        for i in range(len(p1)):
            border1 = max(self.paramSettings.minVals[i],p1[i]-alpha*(p2[i]-p1[i]))
            border2 = max(self.paramSettings.minVals[i],p2[i]+alpha*(p2[i]-p1[i]))
            o.append(np.random.uniform(border1,border2))
        return np.array(o)
    def crossoverSBX(self,p1,p2):
        assert(len(p1)==len(p2))
        beta = np.random.uniform(0,2)
        o1 = []
        o2 = []
        for i in range(len(p1)):
            o1.append(max(self.paramSettings.minVals[i],0.5*((1+beta)*p1[i]+(1-beta)*p2[i])))
            o2.append(max(self.paramSettings.minVals[i],0.5*((1-beta)*p1[i]+(1+beta)*p2[i])))
        return [np.array(o1),np.array(o2)]
    def mutationGaussian(self,p1,mutationRate=1.0,immutableIds=[],width=0.1):
        o = np.copy(p1)
        for i in range(len(p1)):
            if(i not in immutableIds):
                if(mutationRate >= np.random.rand()):
                    o[i] += np.random.normal()*(self.paramSettings.maxVals[i]-self.paramSettings.minVals[i])*width
                    o[i] = max(self.paramSettings.minVals[i],o[i])
        return o
        
class RandomSearch:
    def __init__(self,paramSettings):
        self.paramSettings = paramSettings
    def limitValues(self,p):
        for i in range(len(p)):
            if(p[i] < paramSettings.minVals[i]):
                p[i] = paramSettings.minVals[i]
            if(p[i] > paramSettings.maxVals[i]):
                p[i] = paramSettings.maxVals[i]
        return p
    def optimize(self,startPoint,currentFitness,stepSize,fn,minimumImprovement,maximumSteps=100,dropout=0.0):
        p = np.copy(startPoint)
        d = np.random.standard_normal(size=len(startPoint))
        d = self.paramSettings.denormalizeParameterValues(d)*stepSize*np.random.uniform()
        #Add dropout: Remove part of new direction parameters
        for i in range(len(d)):
            if(dropout > np.random.rand()):
                d[i] = 0
        numOfEvals = 1
        probeFitness = fn(self.limitValues(p+d)) #only positive parameters
        if(probeFitness>currentFitness):
            d = -d
            probeFitness = fn(p+d)
        while(probeFitness<currentFitness and numOfEvals<10):
            if(np.abs(1 - probeFitness/currentFitness) < minimumImprovement):
                break
            currentFitness = probeFitness
            p = self.limitValues(p+d)
            probeFitness = fn(self.limitValues(p+d))
            numOfEvals += 1
        if(probeFitness<currentFitness):
            return [probeFitness,self.limitValues(p+d)]
        else:
            return [currentFitness,p]

class FileLogger:
    def __init__(self,logFileName,mode):
        self.fileName = logFileName
        self.mode = mode
    def write(self,content):
        with open(self.fileName,self.mode) as out:
            out.write(content)

class GenOptCore:
    def __init__(self,paramSettings,conf,fitnessFunction,generations,populationSize,mutationsPerGeneration,stepSize,selectionProbabilities,minimumImprovement,minimumDiversity):
        self.paramSettings = paramSettings
        self.fitnessFunction = fitnessFunction
        self.p = Population(paramSettings,self.fitnessFunction,maxSize=populationSize,mutationsPerGeneration=mutationsPerGeneration)
        self.go = GeneticOperators(paramSettings)
        self.rs = RandomSearch(paramSettings)
        self.generations = generations
        self.stepSize = stepSize
        self.probabilities = selectionProbabilities
        self.minimumImprovement = minimumImprovement
        self.minimumDiversity = minimumDiversity
        
        self.startTime = time.time()
        self.meanFitness = np.array([])
        self.bestFitness = np.array([])
        self.timeLog = [0]
        self.evalLog = [0]
        self.conf = conf
        self.fileLogger = FileLogger(conf.logfilename,"a+")
        self.populationLogger = FileLogger(conf.populationfilename,"w")
        self.config_content = conf.CONFIG_CONTENT
    def geneticOptimization(self,selectionAlgorithm):
        selection = selectionAlgorithm
        newMembers = []
        assert(np.sum(self.probabilities)==1.0)
        probabilities = np.cumsum(self.probabilities)
        for i in range(self.p.mutationsPerGeneration):
            p1 = selection(self.p)
            p2 = selection(self.p)
            rn = np.random.rand()
            if(probabilities[0]>rn):
                offspring = self.go.crossoverBLX(p1.values,p2.values)
                newMembers.append(PopulationMember(offspring,self.fitnessFunction.compute(offspring,"genetic")))
            elif(probabilities[1]>rn):
                offspring1,offspring2 = self.go.crossoverSBX(p1.values,p2.values)
                newMembers.append(PopulationMember(offspring1,self.fitnessFunction.compute(offspring1,"genetic")))
                newMembers.append(PopulationMember(offspring2,self.fitnessFunction.compute(offspring2,"genetic")))
            elif(probabilities[2]>rn):
                offspring = self.go.crossoverUniform(p1.values,p2.values)
                newMembers.append(PopulationMember(offspring,self.fitnessFunction.compute(offspring,"genetic")))
        for i in range(self.p.mutationsPerGeneration):
            p = selection(self.p)
            offspring = self.go.mutationGaussian(p.values)
            newMembers.append(PopulationMember(offspring,self.fitnessFunction.compute(offspring,"genetic")))
        for m in newMembers:
            self.p.addPopulationMember(m)
        self.p.cutPopulation()
    def randomOptimization(self):
        selection = self.go.rankSelection
        for i in range(self.p.mutationsPerGeneration):
            p = selection(self.p)
            offspring = self.rs.optimize(p.values,p.fitness,self.stepSize,self.fitnessFunction.compute,self.minimumImprovement)
            if(p.fitness >  offspring[0]):
                self.p.addPopulationMember(PopulationMember(offspring[1],offspring[0]))
        self.p.cutPopulation()
    def logLineSummary(self,i):
        self.timeLog.append(time.time()-self.startTime)
        self.evalLog.append(self.fitnessFunction.numOfEvaluations)
        return str(str(i+1)+"\t"+
                       "{:.2f}".format(self.bestFitness[-1])+"\t"+
                       "{:.2f}".format(self.meanFitness[-1])+"\t"+
                       str(self.fitnessFunction.numOfEvaluations)+"\t"+
                       str(self.fitnessFunction.numOfGeneticEvaluations)+"\t"+
                       str(self.fitnessFunction.numOfEvaluations-self.fitnessFunction.numOfGeneticEvaluations)+"\t"+
                       "{:.2f}".format(self.timeLog[-1])+"\t"+
                       "{:.2f}".format((self.timeLog[-1])/self.fitnessFunction.numOfEvaluations)+"\t"+
                       "{:.2f}".format(self.bestFitness[-2]-self.bestFitness[-1])+"\t"+
                       "{:.2f}".format((self.bestFitness[-2]-self.bestFitness[-1])/self.fitnessFunction.numOfEvaluations)+"\t"+
                       "{:.2f}".format(self.p.estimateDiversity()))
    def dumpParametersCHARMM(self):
        output = ""
        for i in range(self.paramSettings.dim):
            output += "set " + self.paramSettings.names[i] + " " + "{:.2f}".format(self.p.members[0].values[i]) + "\n"
        return output
    def start(self):
        if(self.p.readPopulation(conf.populationfilename) == False):
            self.p.generateInitialPopulation(self.p.maxSize)
        self.p.sortPopulation()
        self.meanFitness = np.append(self.meanFitness,self.p.getFitnessMean())
        self.bestFitness = np.append(self.bestFitness,self.p.getFitnessBest())
        overhead = "Record started: "+str(datetime.datetime.now())+" (Original fitness from cfg file: "+str(fitnessFunction.compute(self.paramSettings.initVals))+")\n"
        overhead += self.config_content
        overhead +="\nGen\tBestF\tMeanF\tTEvals\tGEvals\tREvals\tt\tt/eval\tDelta\tImprov\tDivers"
        self.fileLogger.write(overhead+"\n")
        print(overhead)
        for i in range(self.generations):
            if(self.p.estimateDiversity()<self.minimumDiversity):
                self.p.cutPopulation(1)
                self.p.generateInitialPopulation(self.p.maxSize)
            self.geneticOptimization(self.go.rankSelection)
            self.meanFitness = np.append(self.meanFitness,self.p.getFitnessMean())
            self.bestFitness = np.append(self.bestFitness,self.p.getFitnessBest())
            if(i != 0 and self.bestFitness[-1] == self.bestFitness[-2]):
                self.randomOptimization()
            logLine = self.logLineSummary(i)
            print(logLine)
            self.fileLogger.write(logLine+"\n")
            self.populationLogger.write("Generation "+str(i+1)+"\n"+self.p.dumpPopulation()+"\n"+self.fitnessFunction.dumpFrequencies(self.p.members[0].values)+"\n"+self.dumpParametersCHARMM())
        print("Gen\tBestF\tMeanF\tTEvals\tGEvals\tREvals\tt\tt/eval\tDelta\tImprov\tDivers\n")
        return self.bestFitness[-1]

conf = Config()
conf.load()

paramSettings = ParameterSettings()
for p in conf.PARAMETER_SETTINGS:
    paramSettings.addParameterSetting(name=p[0],ptype=p[1],pmin=p[2],pmax=p[3],initVal=p[4])
fitnessFunction = FitnessFunction(paramSettings,conf.extern,conf.psf,conf.crd,conf.params,conf.varfile,conf.mdexec,conf.mdinp,conf.mdout,conf.qmout,conf.paramfilename)

core = GenOptCore(paramSettings,conf,fitnessFunction,generations=conf.GENERATIONS,populationSize=conf.POPULATION_SIZE,mutationsPerGeneration=conf.MUTATIONS_PER_GENERATION,stepSize=conf.STEP_SIZE,selectionProbabilities=[conf.CROSSOVER_BLX,conf.CROSSOVER_SBX,conf.CROSSOVER_UNIFORM],minimumImprovement=conf.MINIMUM_IMPROVEMENT,minimumDiversity=conf.MINIMUM_DIVERSITY)

print(core.start())

