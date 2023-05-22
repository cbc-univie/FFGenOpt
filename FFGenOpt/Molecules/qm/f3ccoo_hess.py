import psi4,os,importlib
import numpy as np
import sys

mem='96gb'
cpu=os.cpu_count() # for cpu rs10

#mem='20gb'
#cpu=24 # for cpuongpu rs10

psi4.set_num_threads(cpu)
psi4.set_memory(mem)

outfile = 'f3ccoo_hess.out'
psi4.core.set_output_file(outfile, False)
psi4.core.IOManager.shared_object().set_default_path('/scratch/andras/')

inp_mol = psi4.geometry("""
-1 1
C     -1.60392   -0.46377    0.02162
C     -0.36006    0.41501    0.05794
F     -1.97453   -0.73962   -1.24327
F     -1.40553   -1.63727    0.65260
F     -2.65285    0.13123    0.62146
O      0.71408   -0.17910   -0.39895
O     -0.37609    1.60231    0.32048
""")

mp2_e, mp2_wfn = psi4.optimize('mp2/cc-pvdz', molecule=inp_mol, return_wfn = True)
hess, wfn = psi4.driver.hessian('mp2/cc-pvdz', return_wfn = True)
dipder = wfn.variables().get("CURRENT DIPOLE GRADIENT", None)
if dipder is not None:
    dipder = np.asarray(dipder).T
#print(inp_mol.__dict__)
#print(dir(inp_mol))
hess_arr = np.asarray(hess)
geom = np.asarray(inp_mol.geometry())
masses = np.asarray([inp_mol.mass(i) for i in range(inp_mol.natom())])
ir_labels = inp_mol.irrep_labels()
basis = wfn.basisset()
#wfn.hessian().print_out()

vibinfo, vibtext = psi4.driver.qcdb.vib.harmonic_analysis(hess_arr, geom, masses, basis, ir_labels, dipder)
symbols = [inp_mol.symbol(at) for at in range(inp_mol.natom())]
#print(vibinfo['omega'])
#print(vibinfo['q'])
#print(vibtext)

#psi4.driver.vibanal_wfn(wfn)
with open(outfile, 'a') as f:
    print(psi4.driver.qcdb.vib.print_vibs(vibinfo, normco='q', ncprec=5), file=f)

with open("f3ccoo.molden", 'w') as handle:
    handle.write(psi4.driver.qcdb.vib.print_molden_vibs(vibinfo, symbols, geom, standalone=True))
