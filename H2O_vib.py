from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
#/home/fradion12/ENV/lib64/python3.5/site-packages/pyscf/prop/freq/rhf.py



import numpy
import matplotlib.pyplot as plt

from pyscf.prop import freq
from pyscf import gto, dft

mol = gto.M(atom='H 0.000000  -0.748791  -0.359532; O 0.000000 -0.000000   0.219063; H 0.000000   0.748791  -0.359532',basis='ccpvdz')
#
#Cartesian coordinates (Angstrom)
# Atom        New coordinates             dX        dY        dZ
#   H   0.000000  -0.748791  -0.359532    0.000000  0.000296  0.000040
#   O   0.000000  -0.000000   0.219063    0.000000 -0.000000 -0.000080
#   H   0.000000   0.748791  -0.359532    0.000000 -0.000296  0.000040


#mol.charge = 1
mf = scf.RHF(mol)
#mf.kernel()
mol_eq = optimize(mf)
mol_eq.kernel()


meq = scf.RHF(mol_eq).run()
w, modes = freq.rhf.Freq(meq).kernel() 
print(w)
print(modes)
#atest1,test2 = scf.RHF(mol_eq).Freq(mol_eq)

#mf = dft.RKS(mol_eq).run()

#w, modes = freq.rks.Freq(mf).kernel()
#print(w)
#print(modes)
