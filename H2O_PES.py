from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

#mycc = cc.RCCSD(mf)
#mycc.kernel()
#e,v = mycc.ipccsd(nroots=3)
#print(e)
#print(v)


print(mol.atom_coords())
H1=mol_eq.atom_coords()[0]
O=mol_eq.atom_coords()[1]
H2=mol_eq.atom_coords()[2]
r1 = mol_eq.atom_coords()[1] - mol_eq.atom_coords()[0]
nr1 = numpy.linalg.norm(r1)
print(nr1)
r2 = mol_eq.atom_coords()[1] - mol_eq.atom_coords()[2]
nr2 = numpy.linalg.norm(r2)
print(nr2)
print(H1[0])
print(H1[1])
print(H1[2])

print(O[0])
print(O[1])
print(O[2])

print(H2[0])
print(H2[1])
print(H2[2])

grandR=O[2]-(H1[2])
petitr=numpy.abs(H2[1]-H1[1])
print("grandR="+str(grandR))
print("petitr="+str(petitr))

mol.atom = 'H '+str(H1[0])+' '+str(H1[1])+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2])+'; H '+str(H2[0])+' '+str(H2[1])+' '+str(H2[2])
mol.unit = 'Bohr'
mol.build()
mf = scf.RHF(mol)
mf.kernel()

dR=0.1
dpetitr=0.1
f = open('data_pes.dat','w')
X=[]
dX=[]
ESCF=[]
ERCCSD=[]
EMCSCF=[]
npetitr=13
for idpetitr in range(npetitr):
  ndpetitr=-numpy.floor(npetitr/2)+idpetitr 
  petitr=ndpetitr*dpetitr
  for idR in range(62):
        ndR=idR-10
        print(ndR)
        dOz=ndR*dR
        print(dOz)
        mol.atom = 'H '+str(H1[0])+' '+str(H1[1]-petitr/2)+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2]+dOz)+'; H '+str(H2[0])+' '+str(H2[1]+petitr/2)+' '+str(H2[2])
        mol.build()
        mf = scf.RHF(mol)
        mf.kernel()
#        mycc = cc.CCSD(mf)
#        mycc.kernel()
#        e,v = mycc.ipccsd(nroots=3)
#        f.write(str(dOz)+" "+str(O[2]+dOz-H1[2])+" "+str(mf.energy_tot())+" "+str(mycc.energy())+'\n')
        mc = mcscf.CASCI(mf, 8, 8)
        print('E(CASCI) = %.9g' % mc.casci()[0])
        X.append(O[2]+dOz)
        dX.append(dOz)
        ESCF.append(mf.energy_tot())
        EMCSCF.append(mc.casci()[0])
#        print("mycc.energy()")
#        print(mycc.energy())
#        ERCCSD.append(mf.energy_tot()-mycc.energy())

  #fig = plt.figure()
  plt.plot(X,ESCF,label='SCF, r='+str((H2[1]+petitr/2)-(H1[1]-petitr/2)))
  plt.plot(X,EMCSCF,label='MCSCF, r='+str((H2[1]+petitr/2)-(H1[1]-petitr/2)))
  plt.legend()
  plt.savefig('PES_H2O_r_'+str(((H2[1]+petitr/2)-(H1[1]-petitr/2)))+".png")
  plt.close()
  X=[]
  ESCF=[]
  EMCSCF=[]

plt.legend()
plt.savefig('PES_H2O')
#plt.show()

mol.atom = 'H '+str(H1[0])+' '+str(H1[1])+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2])+'; H '+str(H2[0])+' '+str(H2[1])+' '+str(H2[2])
mol.unit = 'Bohr'
mol.build()
mf = scf.RHF(mol)
mf.kernel()

dR=0.1
X=[]
dX=[]
ESCFion=[]
EMCSCFion=[]
for idR in range(62):
        ndR=idR-10
        print(ndR)
        dOz=ndR*dR
        print(dOz)
        mol.atom = 'H '+str(H1[0])+' '+str(H1[1])+' '+str(H1[2])+'; O '+str(O[0])+' '+str(O[1])+' '+str(O[2]+dOz)+'; H '+str(H2[0])+' '+str(H2[1])+' '+str(H2[2])
        mol.charge = 1
        mol.spin=3
        mol.build()
        mf = scf.UHF(mol)
        mf.kernel()
#        mycc = cc.CCSD(mf)
#        mycc.kernel()
#        e,v = mycc.ipccsd(nroots=3)
#        f.write(str(dOz)+" "+str(O[2]+dOz-H1[2])+" "+str(mf.energy_tot())+" "+str(mycc.energy())+'\n')
        X.append(O[2]+dOz)
        dX.append(dOz)
        ESCFion.append(mf.energy_tot())
        mc = mcscf.CASCI(mf, 8, 7)
        EMCSCFion.append(mc.casci()[0])
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
