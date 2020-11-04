from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
import numpy

import configparser


def getR(mol):
        # center of mass of atom 2 and 3
        R=numpy.linalg.norm(mol.atom_coords()[1]-mol.atom_coords()[0])
        return R


def getpetitr(mol):
        R=numpy.linalg.norm(mol.atom_coords()[2]-mol.atom_coords()[1])
        return R


config = configparser.ConfigParser()
config.read('config.ini')
print(config['DEFAULT']['opt'])
import matplotlib
#['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

matplotlib.use('GTK3Agg')

if config['DEFAULT']['agg']:
   matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

#os.mkdir('CHK')


mol = gto.Mole()
mol.verbose=int(config['DEFAULT']['verbose'])
print(config['DEFAULT']['atom'])
mol.atom = config['DEFAULT']['atom']
mol.basis = config['DEFAULT']['basis']
mol.build()

print('TEST')
atomlist=[]
for ia in range(mol.natm):
        atomlist.append(mol.atom_pure_symbol(ia))

print(atomlist)

mf = scf.RHF(mol)
#mf.kernel()
if config['DEFAULT']['opt']:
        mol_eq = optimize(mf)
	#mol_eq.chkpoint = "CHK/eq.chk"
        mol_eq.kernel()
        print(mol_eq.atom_mass_list())
        # center of mass of atom 2 and 3
        R = getR(mol)
        Req = getR(mol_eq)
        print("R = " + str(R))
        print("R_eq = " + str(Req))
        mol = mol_eq




print(mol_eq)


nval=config['GRID']['nval'].split(',')
minmax=config['GRID']['minmax'].split(',')
dim1vec=numpy.linspace(numpy.float(minmax[0]),numpy.float(minmax[1]),numpy.int(nval[0]))
ESCF=numpy.zeros(numpy.int(nval[0]))





ESCF = numpy.save('SCFPES', ESCF)
print(ESCF)
plt.plot(ESCF)
plt.show()
