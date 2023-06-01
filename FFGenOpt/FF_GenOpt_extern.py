import string
import os
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import parmed

from normalmodeanalysis import NormalModeAnalysis

qmfactor = 0.957

#source: from afmm_read_gaupy.py:
def read_charmm(filename):
    """Takes the name of the CHARMM output file that contains
    normal modes"""

    freq = []
    X = []
    Y = []
    Z = []

    vibmod = "  VIBRATION MODE"
    eigvec = "   EIGENVECTOR:"

    try:
        input = open(filename, 'r')
    except IOError:
        print("Cannot read file " + filename)
        sys.exit(-1)
    x = input.readline()
    while len(x) != 0:
        if x[:len(vibmod)] == vibmod:
            #l = string.split(x)
            l = x.split()
            if len(l[3]) > 12:
                if isinstance(l[3][11:], float):
                    freq.append(float(l[3][11:]))
                else:
                    freq.append(0)
            else:
                freq.append(float(l[4]))
            x = input.readline()
            while x[:len(eigvec)] != eigvec:
                x = input.readline()
            # get next line
            l = input.readline().split()
            while len(l) == 7:
                X.append(float(l[4]))
                Y.append(float(l[5]))
                Z.append(float(l[6]))
                l = input.readline().split()
        x = input.readline()
    input.close()
    
    # number of atoms
    N = int(len(freq) / 3)

    # return results, skipping the zero frequencies
    # this may cause problems with lonepairs as they are also zero and come before the vib freqs
    return freq[6:], X[N * 6:], Y[N * 6:], Z[N * 6:] #original one

#from afmm_read_gaupy.py
def read_gaussian(filename):
    """Takes the name of the QM output file that contains
    normal modes (currently Guassian / Psi4)"""
    freq = []
    X = []
    Y = []
    Z = []
    try:
        input = open(filename, 'r')
    except IOError:
        print("Cannot read file " + filename)
        sys.exit(-1)
    x = input.readline()
    if len(x.split()) > 1 and x.split()[1] == "Gaussian":
        searchstring = " Frequencies"
    elif len(x.split()) == 0:
        searchstring = "  Freq"
    else:
        print("Unknown QM file format. Exiting")
        sys.exit(-1)
    while len(x) != 0:
        if x.find(searchstring) == 0:
            l = x.split()
            for i in range(3):
                freq.append(float(l[i + 2]))
            l = input.readline().split()
            while l[0][0] not in string.digits:
                l = input.readline().split()
            bx1 = []
            by1 = []
            bz1 = []
            bx2 = []
            by2 = []
            bz2 = []
            bx3 = []
            by3 = []
            bz3 = []
            while len(l) != 0:
                bx1.append(float(l[-9]))
                by1.append(float(l[-8]))
                bz1.append(float(l[-7]))
                bx2.append(float(l[-6]))
                by2.append(float(l[-5]))
                bz2.append(float(l[-4]))
                bx3.append(float(l[-3]))
                by3.append(float(l[-2]))
                bz3.append(float(l[-1]))
                l = input.readline().split()
            X = X + bx1 + bx2 + bx3
            Y = Y + by1 + by2 + by3
            Z = Z + bz1 + bz2 + bz3
        x = input.readline()
    input.close()
 
    return freq, X, Y, Z

#source: from afmm_read_gaupy.py:
def Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ):
    """Computes the Merit function"""
    N = int((len(mdfreq)/3) + 2)
    if (N != ((len(qmfreq)/3) + 2)) or (N == 2):
        print("Atoms from MD:", N, ". Atoms from QM:", int((len(qmfreq)/3) + 2))
        print("Different or zero number of atoms. Exiting.")
        sys.exit()
    weight = []
    maxprojidxs = []
    nfreq = len(mdfreq)
    # cost matrix, saves qm and mm dot products, row: MD, col: QM
    costmatrix = np.zeros([nfreq, nfreq]) 
    for mdidx in range(nfreq):
        maxproj = 0.0
        maxprojidx = 0
        for qmidx in range(nfreq):
            # gets the correct values and does the dot product
            mdstart = mdidx * N
            qmstart = qmidx * N
            proj = DotProduct(mdX[mdstart:mdstart+N], mdY[mdstart:mdstart+N], mdZ[mdstart:mdstart+N],
                qmX[qmstart:qmstart+N], qmY[qmstart:qmstart+N], qmZ[qmstart:qmstart+N])
            #build the costmatrix
            try:
                costmatrix[mdidx][qmidx] = proj * min([mdfreq[mdidx]/qmfreq[qmidx],
                    qmfreq[qmidx]/mdfreq[mdidx]])    #1-proj if not maximizue=True
            except ZeroDivisionError:
                print(f"MDFREQ: {mdfreq[mdidx]}, QMFREQ: {qmfreq[mdidx]}")
                costmatrix[mdidx][qmidx] = 0.0
    #hungarian method, to ensure the best 1:1 mapping of QM and MM freqs
    try:
        row_ind, col_ind = linear_sum_assignment(costmatrix, maximize=True)
    except ValueError:
        print("VALUE ERROR IN COSTMATRIX")
        print(costmatrix)
        quit()
    maxprojidxs = col_ind

    for pos, maxprojidx in enumerate(maxprojidxs):
        maxproj = costmatrix[pos][maxprojidx]
        #print(maxproj, max(costmatrix[pos]), (maxproj/max(costmatrix[pos]))*100 ,"%")
        mdidx = pos
        # compute the different weights
        # Attention: "projection" does not work, as maxproj often 0, "frequency" also resuletd in error sometimes
        #if self.weighting == "projection":
        #    if self.TooLow(maxproj):
        #        print("Maximum projection very close to zero. Exiting.")
        #        sys.exit(-1)
        #    weight.append(1/maxproj)
        #elif self.weighting == "frequency":
        #    if self.TooLow(mdfreq[mdidx]):
        #        print("MD frequency very close to zero. Exiting.")
        #        sys.exit(-1)
        #weight.append(1000/mdfreq[mdidx])
        #elif self.weighting == "freqmatch":
        #    md = mdfreq[mdidx]
        #    qm = qmfreq[maxprojidx]*self.qmfactor
        #    weight.append(abs((md - qm)/max(md, qm)))
        #elif self.weighting == "none":
        weight.append(1.0)
    sum = 0.0
    nwsum = 0.0
    mdqmfreq = []
    for i in range(nfreq):
        # compute intermediate result to avoid doing it twice
        t = (mdfreq[i] - qmfactor*qmfreq[maxprojidxs[i]])**2
        # weighted sum of squares
        sum = sum + (weight[i]**2) * t
        # non-weighted sum of squares needed for final output
        nwsum = nwsum + t
        mdqmfreq.append(qmfreq[maxprojidxs[i]])
    return (sum/nfreq)**0.5, (nwsum/nfreq)**0.5, mdfreq, mdqmfreq, maxprojidxs

def DotProduct(mdX, mdY, mdZ, qmX, qmY, qmZ):
        sum = 0.0
        for i in range(len(mdX)):
            sum = sum + mdX[i]*qmX[i] + mdY[i]*qmY[i] + mdZ[i]*qmZ[i]
        return abs(sum)

def RunMD(mdexec, mdinp, mdout):
    """Run MD program, checking for normal termination"""
    cmd = mdexec + " < " + mdinp + " 2> mderror.log 1> " + mdout
    status = os.system(cmd)
    if os.WIFEXITED(status):
        if os.WEXITSTATUS(status) != 0:
            print("MD program returned an error! Not exiting.")
            #print("MD program returned an error! Exiting.")
            #sys.exit(-1)
    else:
        print("MD program signalled! Not exiting.")
        #sys.exit(-1)

def create_context(psffile,crdfile,paramsfile):
    """Create simulation objects in OpenMM"""
    psf = CharmmPsfFile(psffile)
    crd = CharmmCrdFile(crdfile)

    parFiles = ()
    for line in open(paramsfile, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0: parFiles += ( parfile, )

    parmed_prm = parmed.charmm.CharmmParameterSet( *parFiles )
    pmff = parmed.openmm.parameters.OpenMMParameterSet.from_parameterset(parmed_prm,
            unique_atom_types = True)
    pmff.write('parmed_omm.xml', separate_ljforce = True)
    ff = ForceField('parmed_omm.xml')

    params = CharmmParameterSet( *parFiles )
    #params = CharmmParameterSet(paramsfile)

    nonbondedCutoff = 100.*nanometers
    constraintTolerance = 0.000001
    dt = 0.001*picoseconds
    temperature = 300*kelvin
    friction = 1.0/picosecond
    pressure = 1.0*atmospheres
    coll_freq = 10.0/picosecond
    drud_coll_freq = 200.0/picosecond
    drud_temp = 1*kelvin

    platform = Platform.getPlatformByName('CPU')

    topology = []

    positions = crd.positions
    mod = Modeller(psf.topology, positions)
    mod.addExtraParticles(ff)

    topology.append(mod.topology)
    topology.append(psf.topology)

    system = []
    system.append(ff.createSystem(mod.topology, nonbondedCutoff=nonbondedCutoff))
    system.append(psf.createSystem(params, nonbondedCutoff=nonbondedCutoff))

    integrator = []
    integrator.append(DrudeNoseHooverIntegrator(temperature, coll_freq, drud_temp, drud_coll_freq, dt))
    integrator[0].setConstraintTolerance(constraintTolerance)
    integrator.append(NoseHooverIntegrator(temperature, coll_freq, dt))
    integrator[1].setConstraintTolerance(constraintTolerance)
    context = []
    context.append(Context(system[0], integrator[0], platform))
    context[0].setPositions(mod.positions)
    context.append(Context(system[1], integrator[1], platform))
    context[1].setPositions(crd.positions)

    return context, mod, system, integrator, psf, topology

def get_varnames(streamfile):
    """Get new force constants from stream file"""
    varnames = {}
    parfile = open(streamfile, 'r')
    x = parfile.readline()
    while len(x) != 0:
        varnames[x.split()[1]] = float(x.split()[2])
        x = parfile.readline()
    parfile.close()

    return varnames

def to_change(vars, varfile, psf, streamfile, varnames):
    """Determine which parameters are being optimized"""
    if varnames is None:
        varnames = get_varnames(streamfile)

    def explicit_params(pnames, param_names):
        """Substitute wildcards by explicit atomtypes"""
        if len(pnames) != 6:
            print("Currently only dihedral/improper wildcards supported. Exiting.")
            sys.exit(-1)
        expl_attypes = []
        if pnames[1] == 'X' and pnames[4] == 'X':
            for i in range(len(psf.dihedral_list)):
                if psf.dihedral_list[i].atom2.attype == pnames[2]\
                    and psf.dihedral_list[i].atom3.attype == pnames[3]:
                    expl_attypes.append((pnames[0], psf.dihedral_list[i].atom1.attype,
                    psf.dihedral_list[i].atom2.attype, psf.dihedral_list[i].atom3.attype,
                    psf.dihedral_list[i].atom4.attype, pnames[5]))

                elif psf.dihedral_list[i].atom2.attype == pnames[3]\
                    and psf.dihedral_list[i].atom3.attype == pnames[2]:
                    expl_attypes.append((pnames[0], psf.dihedral_list[i].atom4.attype,
                    psf.dihedral_list[i].atom3.attype, psf.dihedral_list[i].atom2.attype,
                    psf.dihedral_list[i].atom1.attype, pnames[5]))


        if pnames[2] == 'X' and pnames[3] == 'X':
            for i in range(len(psf.improper_list)):
                if psf.improper_list[i].atom1.attype == pnames[1]\
                    and psf.improper_list[i].atom4.attype == pnames[4]:
                    expl_attypes.append((pnames[0], psf.improper_list[i].atom1.attype,
                    psf.improper_list[i].atom2.attype, psf.improper_list[i].atom3.attype,
                    psf.improper_list[i].atom4.attype, pnames[5]))

                elif psf.improper_list[i].atom1.attype == pnames[4]\
                    and psf.improper_list[i].atom4.attype == pnames[1]:
                    expl_attypes.append((pnames[0], psf.improper_list[i].atom4.attype,
                    psf.improper_list[i].atom3.attype, psf.improper_list[i].atom2.attype,
                    psf.improper_list[i].atom1.attype, pnames[5]))

        clean_attypes = list(set(expl_attypes))
        list_attypes  = [list(a) for a in clean_attypes]
        for i in range(len(list_attypes)):
            list_attypes[i][0] += f'_{i}'
            varnames[list_attypes[i][0]] = varnames[clean_attypes[i][0]]
            param_names.append(list_attypes[i])

    param_names = []
    bond_types = {}
    angle_types = {}
    dih_types = {}
    param_orig = open(varfile, 'r')
    x = param_orig.readline()
    while len(x) != 0:
        if len(x.split()) > 0 and x.split()[0].startswith("!"):
            x = param_orig.readline()
            continue
        x = x.split("!", 1)[0]
        for i in x.split():
            if i[1:] in vars:
                last_elem = x.split().index(i)
                pnames = [i[1:]]
                for j in range(last_elem):
                    pnames.append(x.split()[j])
                if last_elem == 4:
                    pnames.append(x.split()[5])
                if 'X' in pnames:
                    explicit_params(pnames, param_names)
                param_names.append(pnames)
        x = param_orig.readline()
    #print(param_names)
    param_orig.close()
    
    for i in range(len(param_names)):
        if len(param_names[i]) == 3:
            bond_types[param_names[i][0]] = param_names[i][1:3]
        elif len(param_names[i]) == 4:
            angle_types[param_names[i][0]] = param_names[i][1:4]
        elif len(param_names[i]) == 6:
            dih_types[param_names[i][0]] = param_names[i][1:6]
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
                and psf.bond_list[j].atom2.attype == bond_types[i][1])\
                or (psf.bond_list[j].atom1.attype == bond_types[i][1]\
                and psf.bond_list[j].atom2.attype == bond_types[i][0]):
    
                bonds_to_change[(psf.bond_list[j].atom1.idx, psf.bond_list[j].atom2.idx)] = i
 
    if hasattr(psf, "urey_bradley_list"):
        for i in bond_types:
            for j in range(len(psf.urey_bradley_list)):
                if (psf.urey_bradley_list[j].atom1.attype == bond_types[i][0]
                    and psf.urey_bradley_list[j].atom2.attype == bond_types[i][1])\
                    or (psf.urey_bradley_list[j].atom1.attype == bond_types[i][1]\
                    and psf.urey_bradley_list[j].atom2.attype == bond_types[i][0]):
        
                    bonds_to_change[(psf.bond_list[j].atom1.idx, psf.bond_list[j].atom2.idx)] = i
        #print(bonds_to_change)
    
    angles_to_change = {}
    for i in angle_types:
        for j in range(len(psf.angle_list)):
            if (psf.angle_list[j].atom1.attype == angle_types[i][0]
                and psf.angle_list[j].atom2.attype == angle_types[i][1]\
                and psf.angle_list[j].atom3.attype == angle_types[i][2])\
                or (psf.angle_list[j].atom1.attype == angle_types[i][2]\
                and psf.angle_list[j].atom2.attype == angle_types[i][1]\
                and psf.angle_list[j].atom3.attype == angle_types[i][0]):
    
                angles_to_change[(psf.angle_list[j].atom1.idx,
                    psf.angle_list[j].atom2.idx, psf.angle_list[j].atom3.idx)] = i
    #print(angles_to_change)

    dihs_to_change = {}
    for i in dih_types:
        for j in range(len(psf.dihedral_list)):
            if (psf.dihedral_list[j].atom1.attype == dih_types[i][0]
                and psf.dihedral_list[j].atom2.attype == dih_types[i][1]\
                and psf.dihedral_list[j].atom3.attype == dih_types[i][2]\
                and psf.dihedral_list[j].atom4.attype == dih_types[i][3])\
                or (psf.dihedral_list[j].atom1.attype == dih_types[i][3]\
                and psf.dihedral_list[j].atom2.attype == dih_types[i][2]\
                and psf.dihedral_list[j].atom3.attype == dih_types[i][1]\
                and psf.dihedral_list[j].atom4.attype == dih_types[i][0]):
    
                dihs_to_change[(psf.dihedral_list[j].atom1.idx, psf.dihedral_list[j].atom2.idx,
                    psf.dihedral_list[j].atom3.idx, psf.dihedral_list[j].atom4.idx,
                    int(dih_types[i][4]))] = i
    
    for i in dih_types:
        for j in range(len(psf.improper_list)):
            if (psf.improper_list[j].atom1.attype == dih_types[i][0]
                and psf.improper_list[j].atom2.attype == dih_types[i][1]\
                and psf.improper_list[j].atom3.attype == dih_types[i][2]\
                and psf.improper_list[j].atom4.attype == dih_types[i][3])\
                or (psf.improper_list[j].atom1.attype == dih_types[i][3]\
                and psf.improper_list[j].atom2.attype == dih_types[i][2]\
                and psf.improper_list[j].atom2.attype == dih_types[i][1]\
                and psf.improper_list[j].atom4.attype == dih_types[i][0]):
    
                dihs_to_change[(psf.improper_list[j].atom1.idx, psf.improper_list[j].atom2.idx,
                    psf.improper_list[j].atom3.idx, psf.improper_list[j].atom4.idx)] = i
    
    #print(dihs_to_change)
    return bonds_to_change, angles_to_change, dihs_to_change, varnames

def update_context(system, context, varnames, bonds_to_change, angles_to_change, dihs_to_change):
    """Update force constants in OpenMM context"""
    for f in system.getForces():
        if type(f).__name__ == "HarmonicBondForce":
            for bidx in range(f.getNumBonds()):
                bond = f.getBondParameters(bidx)
                bondinfo = (bond[0], bond[1])
                l  = bond[2]
                if bondinfo in bonds_to_change:
                    f.setBondParameters(bidx, bondinfo[0], bondinfo[1],
                        l, varnames[bonds_to_change[bondinfo]]*kilocalories*angstrom**-2*mole**-1)
                    f.updateParametersInContext(context)
    
        elif type(f).__name__ == "HarmonicAngleForce":
            for aidx in range(f.getNumAngles()):
                angle = f.getAngleParameters(aidx)
                angleinfo = (angle[0], angle[1], angle[2])
                l = angle[3]
                if angleinfo in angles_to_change:
                    f.setAngleParameters(aidx, angleinfo[0], angleinfo[1], angleinfo[2],
                        l, varnames[angles_to_change[angleinfo]]*kilocalories*mole**-1*radian**-2)
                    f.updateParametersInContext(context)
    
        elif type(f).__name__ == "PeriodicTorsionForce":
            for didx in range(f.getNumTorsions()):
                dih = f.getTorsionParameters(didx)
                dihinfo = (dih[0], dih[1], dih[2], dih[3], dih[4])
                delta = dih[5]
                if dihinfo in dihs_to_change:
                    f.setTorsionParameters(didx, dihinfo[0], dihinfo[1], dihinfo[2],
                        dihinfo[3], dihinfo[4], delta,
                        varnames[dihs_to_change[dihinfo]]*kilocalories*mole**-1)
                    f.updateParametersInContext(context)
    
        elif type(f).__name__ == "CustomTorsionForce":
            for iidx in range(f.getNumTorsions()):
                imp = f.getTorsionParameters(iidx)
                impinfo = (imp[0], imp[1], imp[2], imp[3])
                theta0 = imp[4][1]
                if impinfo in dihs_to_change:
                    f.setTorsionParameters(iidx, impinfo[0], impinfo[1],
                        impinfo[2], impinfo[3],
                        (varnames[dihs_to_change[impinfo]], theta0))
                    f.updateParametersInContext(context)

def normal_mode(system, integrator, context, topology):
    """Compute frequencies and normal modes"""
    nma = NormalModeAnalysis(topology[0], system[0],
            integrator[0],
            context[0].getState(getPositions=True).getPositions(asNumpy=True),
            CPUOnly=False)
    nma.CPUPreMinimization()
    #nma.CPUMinimizationCycle()

    heavy_atoms = []
    for i,j in enumerate(topology[0].atoms()):
        if j.element is not None:
            heavy_atoms.append(i)

    minState = nma.CPUSimulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
    minPos   = minState.getPositions(asNumpy=True)[heavy_atoms]
    context[1].setPositions(minPos)
    nma2 = NormalModeAnalysis(topology[1], system[1],
            integrator[1],
            context[1].getState(getPositions=True).getPositions(asNumpy=True),
            CPUOnly=True)

    nma2.CalculateNormalModes()
    vib_spec = []
    cycles = 0
    while float(nma2.VibrationalSpectrum[0]._value) != float(nma2.VibrationalSpectrum[0]._value) and cycles < 1000:
        nma2.CPUMinimizationCycle()
        nma2.CalculateNormalModes()
        cycles += 1
    for i in nma2.VibrationalSpectrum:
        vib_spec.append(float(i._value))
    return vib_spec, nma2.NormalModes[6:,0::3].flatten(), nma2.NormalModes[6:,1::3].flatten(),nma2.NormalModes[6:,2::3].flatten()
