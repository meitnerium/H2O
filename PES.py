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


mol = gto.Mole()
mol.verbose=int(config['DEFAULT']['verbose'])
print(config['DEFAULT']['atom'])
mol.atom = config['DEFAULT']['atom']
mol.basis = config['DEFAULT']['basis']
mol.build()

mf = scf.RHF(mol)
#mf.kernel()
if config['DEFAULT']['opt']:
        mol_eq = optimize(mf)
        mol_eq.kernel()
        print(mol_eq.atom_mass_list())
        # center of mass of atom 2 and 3
        COMA23 = mol_eq.atom_coords()[1]*mol_eq.atom_mass_list()[1]-mol_eq.atom_coords()[2]*mol_eq.atom_mass_list()[2]
        R=numpy.linalg.norm(mol_eq.atom_coords()[0]-COMA23)
        print("R = " + str(R))
        R = getR(mol)
        Req = getR(mol_eq)
        print("R = " + str(R))
        print("R_eq = " + str(Req))
        petitr = getpetitr(mol)
        petitreq = getpetitr(mol_eq)
        print("r = " + str(petitr))
        print("r_eq = " + str(petitreq))
        mol = mol_eq


print(mol_eq)

#mycc = cc.RCCSD(mf)
#mycc.kernel()
#e,v = mycc.ipccsd(nroots=3)
#print(e)
#print(v)

nval=config['GRID']['nval'].split(',')
minmax=config['GRID']['minmax'].split(',')
Rvec=numpy.linspace(numpy.float(minmax[0]),numpy.float(minmax[1]),numpy.int(nval[0]))
petitrvec=numpy.linspace(numpy.float(minmax[2]),numpy.float(minmax[3]),numpy.int(nval[1]))
print(minmax[0])

print(len(Rvec))
print(len(petitrvec))
ESCF=numpy.zeros([len(Rvec),len(petitrvec)])
ndim1=0
molsurf=list(nval)

for dim1 in dim1vec:
	ndim2=0
	''' for Jabcobi variables, dim1 is R, the distance between atom 1 and the centre of mass of atom 2 and 3 '''
	A1=[0,0,0]
	A2=[0,0,-dim1]
	A3=[0,0,-dim1]
	ndimm2=0
	for dim2 in dim2vec:
		''' for Jabcobi variables, dim2 is R, the distance between atom 2 and 3 '''
		A2[1]=-dim2/2
		A3[1]=dim2/2
		mol.atom = mol.atom[0][0]+' '+str(A1[0])+' '+str(A1[2])+' '+str(A1[3])+' ; '+\
                   mol.atom[1][0]+'H '+str(A1[0])+' '+str(A1[2])+' '+str(A1[3])+' ; '+\
                   mol.atom[2][0]+'H '+str(A1[0])+' '+str(A1[2])+' '+str(A1[3])+' '
		if dim1 == 0 and dim2 == 0:
			print("first calculation")
            	elif npetitr == 0:
               		import shutil
               		oldchk="CHK/"+str(nr-1)+"_"+str(npetitr)+".chk"
               		shutil.copy2(oldchk,"CHK/"+str(nr)+"_"+str(npetitr)+".chk")
            	else:
               		import shutil
               		oldchk="CHK/"+str(nr)+"_"+str(npetitr-1)+".chk"
               		shutil.copy2(oldchk,"CHK/"+str(nr)+"_"+str(npetitr)+".chk")

		mol.buid()
		molsurf[dim1,dim2]=mol
		mf = scf.RHF(mol)
            	mf.chkfile="CHK/"+str(ndim1)+"_"+str(ndim2)+".chk"
            	if nr == 0 and npetitr == 0:
               		print("first calculation")
            	else:
               		mf.init_guess = 'chkfile'
            mf.kernel()
            ESCF[nr,npetitr] = mf.energy_tot()
            npetitr = npetitr +1
        #plt.plot(petitrvec,ESCF[nr,:])
        #plt.show()
        #plt.close()
        nr = nr + 1

#        f.write(str(dOz)+" "+str(O[2]+dOz-H1[2])+" "+str(mf.energy_tot())+" "+str(mycc.energy())+'\n')
#        mc = mcscf.CASCI(mf, 8, 8)
        #print('E(CASCI) = %.9g' % mc.casci()[0])
#        ESCF.append(mf.energy_tot())
#        EMCSCF.append(mc.casci()[0])
#        print("mycc.energy()")
#        print(mycc.energy())
#        ERCCSD.append(mf.energy_tot()-mycc.energy())

  #fig = plt.figure()
#  plt.plot(X,ESCF,label='SCF, r='+str((H2[1]+petitr/2)-(H1[1]-petitr/2)))
#  plt.plot(X,EMCSCF,label='MCSCF, r='+str((H2[1]+petitr/2)-(H1[1]-petitr/2)))
#  plt.legend()
#  plt.savefig('PES_H2O_r_'+str(((H2[1]+petitr/2)-(H1[1]-petitr/2)))+".png")
#  plt.close()
#  X=[]
#  ESCF=[]
#  EMCSCF=[]

#plt.legend()
#plt.savefig('PES_H2O')
#plt.show()

#mol.atom = 'H '+str(H1[0])+' '+str(H1[1])+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2])+'; H '+str(H2[0])+' '+str(H2[1])+' '+str(H2[2])
#mol.unit = 'Bohr'
#mol.build()
#mf = scf.RHF(mol)
#mf.kernel()


