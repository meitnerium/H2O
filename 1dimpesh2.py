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
#['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

import matplotlib
#matplotlib.use('WX')
import matplotlib.pyplot as plt
import os

#os.mkdir('CHK')


mol = gto.Mole()
mol.verbose=int(config['DEFAULT']['verbose'])
print(config['DEFAULT']['atom'])
mol.atom = config['DEFAULT']['atom']
mol.basis = config['DEFAULT']['basis']
mol.charge = int(config['DEFAULT']['charge'])
mol.spin = int(config['DEFAULT']['spin'])
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





ndim1=0
for dim1 in dim1vec:
	A1=[0,0,-dim1/2]
	A2=[0,0,-dim1/2]
	mol.atom = atomlist[0]+' '+str(A1[0])+' '+str(A1[2])+' '+str(-dim1/2)+' ; '+\
        atomlist[1]+' '+str(A1[0])+' '+str(A1[2])+' '+str(dim1/2)
	newchk="CHK/"+str(ndim1)+".chk"
	print(mol.atom)
	if ndim1 == 0:
		print("first calculation")
	else:
		print("not first calculation")
		#import shutil
		#oldchk="CHK/"+str(ndim1-1)+".chk"
		#shutil.copy2(oldchk,newchk)
	mol.build()
	mf = scf.RHF(mol)
	mf.chkfile=newchk
	if dim1 == 0:
		print("first calculation")
	else:
      mf.init_guess = 'chkfile'
   
   mf.kernel()
   
   print(mf.mo_energy)
   
   exit(0)
   
   ESCF[ndim1] = mf.energy_tot()
	ndim1 = ndim1 +1

print(ESCF)
numpy.save('SCFPES', ESCF)
plt.plot(dim1vec,ESCF)
minval=numpy.min(ESCF)
maxval=ESCF[-1]
plt.axis([dim1vec[0],dim1vec[-1],minval-(0.1*(maxval-minval)),maxval+0.25*(maxval-minval)])


plt.savefig('test.png')
plt.show()
plt.close()
