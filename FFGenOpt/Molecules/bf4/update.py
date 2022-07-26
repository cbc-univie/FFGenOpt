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
constraints = False
rigidWater = False
constraintTolerance = 0.000001

# Integration Options

dt = 0.001*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres

# Simulation Options

platform = Platform.getPlatformByName('CPU')
#platformProperties = {'Precision': 'mixed'}

topology = psf.topology
positions = crd.positions
system = psf.createSystem(params, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater)
integrator = NoseHooverIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

varnames = {}
parfile = open("../../FF_GenOpt_results/bf4/parameters.str", 'r')
x = parfile.readline()
while len(x) != 0:
    varnames[x.split()[1]] = float(x.split()[2])
    x = parfile.readline()
#print(varnames)
parfile.close()

param_names = []
bond_types = {}
angle_types = {}
dih_types = {}
param_orig = open("/site/raid4/student6/johnny/sandbox/clap_aug2016_pol.prm")
x = param_orig.readline()
while len(x) != 0:
    if len(x.split()) > 0 and x.split()[0].startswith("!"):
        x = param_orig.readline()
        continue
    x = x.split("!", 1)[0]
    for i in x.split():
        if i[1:] in [a[0] for a in varnames]:
            last_elem = x.split().index(i)
            pnames = [i[1:]]
            for j in range(last_elem):
                pnames.append(x.split()[j])
            param_names.append(pnames)
    x = param_orig.readline()
#print(param_names)
param_orig.close()

for i in range(len(param_names)):
    if len(param_names[i]) == 3:
        bond_types[param_names[i][0]] = param_names[i][1:3]
    elif len(param_names[i]) == 4:
        angle_types[param_names[i][0]] = param_names[i][1:4]
    elif len(param_names[i]) == 5:
        dih_types[param_names[i][0]] = param_names[i][1:5]
    else:
        print("Paramter", i, "not recognized as bond, angle, or dihedral. Exiting")
        sys.exit(-1)
        
#print(bond_types)
#print(angle_types)
#print(dih_types)

bonds_to_change = {}
for i in bond_types:
    for j in range(len(psf.bond_list)):
        if (psf.bond_list[j].atom1.attype == bond_types[i][0]
            and psf.bond_list[j].atom2.attype == bond_types[i][1])
            or (psf.bond_list[j].atom1.attype == bond_types[i][1]
            and psf.bond_list[j].atom2.attype == bond_types[i][0]):

            bonds_to_change[(j, psf.bond_list[j].atom1.idx, psf.bond_list[j].atom2.idx)] = i

        elif (psf.urey_bradley_list[j].atom1.attype == bond_types[i][0]
            and psf.urey_bradley_list[j].atom2.attype == bond_types[i][1])
            or (psf.urey_bradley_list[j].atom1.attype == bond_types[i][1]
            and psf.urey_bradley_list[j].atom2.attype == bond_types[i][0]):

            bonds_to_change[(j, psf.bond_list[j].atom1.idx, psf.bond_list[j].atom2.idx)] = i
#print(bonds_to_change)

angles_to_change = {}
for i in angle_types:
    for j in range(len(psf.angle_list)):
        if (psf.angle_list[j].atom1.attype == angle_types[i][0]
            and psf.angle_list[j].atom2.attype == angle_types[i][1]
            and psf.angle_list[j].atom3.attype == angle_types[i][2])
            or (psf.angle_list[j].atom1.attype == angle_types[i][2]
            and psf.angle_list[j].atom2.attype == angle_types[i][1]
            and psf.angle_list[j].atom3.attype == angle_types[i][0]):

            angles_to_change[(j, psf.angle_list[j].atom1.idx,
                psf.angle_list[j].atom2.idx, psf.angle_list[j].atom3.idx)] = i
#print(angles_to_change)

dihs_to_change = {}
for i in dih_types:
    for j in range(len(psf.dihedral_list)):
        if (psf.dihedral_list[j].atom1.attype == dih_types[i][0]
            and psf.dihedral_list[j].atom2.attype == dih_types[i][1]
            and psf.dihedral_list[j].atom3.attype == dih_types[i][2]
            and psf.dihedral_list[j].atom3.attype == dih_types[i][3])
            or (psf.dihedral_list[j].atom1.attype == dih_types[i][3]
            and psf.dihedral_list[j].atom2.attype == dih_types[i][2]
            and psf.dihedral_list[j].atom2.attype == dih_types[i][1]
            and psf.dihedral_list[j].atom3.attype == dih_types[i][0]):

            dihs_to_change[(j, psf.dihedral_list[j].atom1.idx, psf.dihedral_list[j].atom2.idx,
                psf.dihedral_list[j].atom3.idx, psf.dihedral_list[j].atom4.idx)] = i

        elif (psf.improper_list[j].atom1.attype == dih_types[i][0]
            and psf.improper_list[j].atom2.attype == dih_types[i][1]
            and psf.improper_list[j].atom3.attype == dih_types[i][2]
            and psf.improper_list[j].atom3.attype == dih_types[i][3])
            or (psf.improper_list[j].atom1.attype == dih_types[i][3]
            and psf.improper_list[j].atom2.attype == dih_types[i][2]
            and psf.improper_list[j].atom2.attype == dih_types[i][1]
            and psf.improper_list[j].atom3.attype == dih_types[i][0]):

            dihs_to_change[(j, psf.improper_list[j].atom1.idx, psf.improper_list[j].atom2.idx,
                psf.improper_list[j].atom3.idx, psf.improper_list[j].atom4.idx)] = i

#print(dihs_to_change)

for f in system.getForces():
    if type(f).__name__ == "HarmonicBondForce":
        for bidx in range(f.getNumBonds()):
            bond = f.getBondParameters(bidx)
            bondinfo = (bidx, bond[0], bond[1])
            l  = bond[2]
            if bondinfo in bonds_to_change:
                f.setBondParameters(bondinfo[0], bondinfo[1], bondinfo[2],
                    l, varnames[bonds_to_change[bondinfo]]*kilocalories*angstrom**-2*mole**-1)
                f.updateParametersInContext(simulation.context)

    elif type(f).__name__ == "HarmonicAngleForce":
        for aidx in range(f.getNumAngles()):
            angle = f.getAngleParameters(aidx)
            angleinfo = (aidx, angle[0], angle[1], angle[2])
            l = angle[3]
            if angleinfo in angles_to_change:
                f.setAngleParameters(angleinfo[0], angleinfo[1], angleinfo[2], angleinfo[3],
                    l, varnames[angles_to_change[angleinfo]]*kilocalories*mole**-1*radian**-2)
                f.updateParametersInContext(simulation.context)

    elif type(f).__name__ == "PeriodicTorsionForce":
        for didx in range(f.getNumTorsions()):
            dih = f.getTorsionParameters(didx)
            dihinfo = (didx, dih[0], dih[1], dih[2], dih[3])
            n = dih[4]
            delta = dih[5]
            if dihinfo in dihs_to_change:
                f.setTorsionParameters(dihinfo[0], dihinfo[1], dihinfo[2], dihinfo[3],
                    dihinfo[4], n, delta,
                    varnames[dihs_to_change[dihinfo]]*kilocalories*mole**-1)
                f.updateParametersInContext(simulation.context)

    elif type(f).__name__ == "CustomTorsionForce":
        for iidx in range(f.getNumTorsions()):
            imp = f.getTorsionParameters(iidx)
            impinfo = (iidx, imp[0], imp[1], imp[2], imp[3])
            theta0 = imp[4][2]
            if impinfo in dihs_to_change:
                f.setTorsionParameters(impinfo[0], impinfo[1], impinfo[2],
                    impinfo[3], impinfo[4],
                    (varnames[dihs_to_change[dihinfo]]*kilocalories*mole**-1, theta0))
                f.updateParametersInContext(simulation.context)
