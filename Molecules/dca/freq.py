# This script was generated by OpenMM-Setup on 2021-07-22.

import sys
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

# Input Files

psf = CharmmPsfFile('bf4.psf')
crd = CharmmCrdFile('bf4_new.cor')
params = CharmmParameterSet('clap_aug2016_pol.rtf', 'clap_aug2016_pol.prm')

# System Configuration

nonbondedCutoff = 100.*nanometers
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001

# Integration Options

dt = 0.001*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres

# Simulation Options

platform = Platform.getPlatformByName('CPU')
#platformProperties = {'Precision': 'mixed'}

#xtl = 58.0*angstroms
#psf.setBox(xtl,xtl,xtl)
topology = psf.topology
positions = crd.positions
system = psf.createSystem(params, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater)
integrator = NoseHooverIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

steps = 100000
equilibrationSteps = 500
#checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)


# Prepare the Simulation

print('Building system...')
simulation.minimizeEnergy()

nma = NormalModeAnalysis(topology, system, integrator, positions, CPUOnly=True)
nma.CalculateNormalModes()
print(nma.VibrationalSpectrum)
print(nma.TransRotFreq)
print(nma.NormalModes[6:])
#colors = ["black",]*9
#nma.PlotVibrationalSpectrum(colorStr=colors)
