from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *


psf = CharmmPsfFile("Molecules/bf4/bf4.psf")
crd = CharmmCrdFile("Molecules/bf4/bf4_new.cor")
#params = CharmmParameterSet("/home/andras/git_test/ir_optimization/Molecules/bf4/clap_aug2016_pol.rtf", "/home/andras/git_test/ir_optimization/Molecules/bf4/clap_aug2016_pol.prm")
params = CharmmParameterSet("Molecules/bf4/toppar.str")

nonbondedCutoff = 100.*nanometers
constraintTolerance = 0.000001
dt = 0.001*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
platform = Platform.getPlatformByName('CPU')

topology = psf.topology
positions = crd.positions
system = psf.createSystem(params, nonbondedCutoff=nonbondedCutoff)

integrator = NoseHooverIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

