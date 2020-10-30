from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
import numpy
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

mol = gto.Mole()
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
nr=0
for grandr in Rvec:
	npetitr=0
	for petitr in petitrvec:
		A1Z = grandr*mol_eq.atom_mass_list()[0]/(mol_eq.atom_mass_list()[0]+mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2])
		A23Z = grandr*(mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2])/(2*(mol_eq.atom_mass_list()[0]+mol_eq.atom_mass_list()[1]+mol_eq.atom_mass_list()[2]))
		A23Y = petitr/2
		mol.atom = 'O 0.0 0.0 '+str(A1Z)+' ; H 0.0 -'+str(A23Y)+' '+str(-A23Z)+' ; H 0.0 '+str(A23Y)+' '+str(-A23Z)
		mol.unit = 'Bohr'
		mol.build()
		mf = scf.RHF(mol)
		mf.kernel()
		ESCF[nr,npetitr] = mf.energy_tot() 
		npetitr = npetitr +1
	plt.plot(petitrvec,ESCF[nr,:])
	plt.show()
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

#dR=0.1
#X=[]
#dX=[]
#ESCFion=[]
#EMCSCFion=[]
#for idR in range(62):
#        ndR=idR-10
#        print(ndR)
#        dOz=ndR*dR
#        print(dOz)
#        mol.atom = 'H '+str(H1[0])+' '+str(H1[1])+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2]+dOz)+'; H '+str(H2[0])+' '+str(H2[1])+' '+str(H2[2])
#        mol.charge = 1
#        mol.spin=3
#        mol.build()
#        mf = scf.UHF(mol)
#        mf.kernel()
##        mycc = cc.CCSD(mf)
#        mycc.kernel()
#        e,v = mycc.ipccsd(nroots=3)
#        f.write(str(dOz)+" "+str(O[2]+dOz-H1[2])+" "+str(mf.energy_tot())+" "+str(mycc.energy())+'\n')
#        X.append(O[2]+dOz)
#        dX.append(dOz)
#        ESCFion.append(mf.energy_tot())
#        mc = mcscf.CASCI(mf, 8, 7)
#        EMCSCFion.append(mc.casci()[0])
#        print("mycc.energy()")
#        print(mycc.energy())
#        ERCCSD.append(mf.energy_tot()-mycc.energy())

#plt.plot(dX,ESCF,label=r'SCF $H_2O$')
#plt.plot(dX,ESCFion,label=r'SCF $H_2O^+$')
#plt.plot(dX,EMCSCF,label=r'MCSCF $H_2O$')
#plt.plot(dX,EMCSCFion,label=r'MCSCF $H_2O^+$')
#plt.xlabel('R (bohr)')
#plt.ylabel('Energy (Hartree)')
#plt.plot(dX,ERCCSD)
#plt.legend()
#plt.savefig('PES_H2O')
#plt.show()

numpy.save('SCFPES', ESCF)
