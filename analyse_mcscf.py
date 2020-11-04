#MCSCF0_0.chk
from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
import numpy

import configparser
def getR(mol):
        # center of mass of atom 2 and 3
        COMA23 = mol.atom_coords()[1]*mol.atom_mass_list()[1]-mol.atom_coords()[2]*mol.atom_mass_list()[2]
        R=numpy.linalg.norm(mol.atom_coords()[0]-COMA23)
        return R


def getpetitr(mol):
        R=numpy.linalg.norm(mol.atom_coords()[2]-mol.atom_coords()[1])
        return R


config = configparser.ConfigParser()
config.read('config.ini')
print(config['DEFAULT']['opt'])
if config['DEFAULT']['agg']:
   import matplotlib
   matplotlib.use('Agg')
import matplotlib.pyplot as plt





nval=config['GRID']['nval'].split(',')
minmax=config['GRID']['minmax'].split(',')
Rvec=numpy.linspace(numpy.float(minmax[0]),numpy.float(minmax[1]),numpy.int(nval[0]))
petitrvec=numpy.linspace(numpy.float(minmax[2]),numpy.float(minmax[3]),numpy.int(nval[1]))
print(minmax[0])

print(len(Rvec))
print(len(petitrvec))
EMCSCF=numpy.zeros([len(Rvec),len(petitrvec)])
nr=0












mol = gto.Mole()
mol.verbose=int(config['DEFAULT']['verbose'])
print(config['DEFAULT']['atom'])
mol.atom = config['DEFAULT']['atom']
mol.basis = config['DEFAULT']['basis']
mol.build()
mol_eq = mol


grandr=Rvec[0]
npetitr=0
for petitr in petitrvec:
            A1Z = grandr*mol_eq.atom_mass_list()[0]/(mol_eq.atom_mass_list()[0]+mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2])
            A23Z = grandr*(mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2])/(2*(mol_eq.atom_mass_list()[0]+mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2]))
            A23Y = petitr/2
            mol.atom = 'O 0.0 0.0 '+str(A1Z)+' ; H 0.0 -'+str(A23Y)+' '+str(-A23Z)+' ; H 0.0 '+str(A23Y)+' '+str(-A23Z)
            mol.unit = 'Bohr'
            import shutil
            newchk="CHK/MCSCF"+str(nr)+"_"+str(npetitr)+".chk"

            mol.build()
            mf = scf.RHF(mol)
            mf.chkfile=newchk
            mf.init_guess = 'chkfile'
            mf.kernel()
            mc = mcscf.CASCI(mf, 16, 8)
            npetitr = npetitr +1
            EMCSCF[nr,npetitr] = mc.casci()[0]

numpy.save('MCSCFPES', EMCSCF)

