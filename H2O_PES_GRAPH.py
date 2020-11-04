from pyscf import gto,scf, ao2mo, cc, mcscf
from pyscf.geomopt.berny_solver import optimize
import numpy
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import configparser


config = configparser.ConfigParser()
config.read('config.ini')

nval=config['GRID']['nval'].split(',')
minmax=config['GRID']['minmax'].split(',')
Rvec=numpy.linspace(numpy.float(minmax[0]),numpy.float(minmax[1]),numpy.int(nval[0]))
petitrvec=numpy.linspace(numpy.float(minmax[2]),numpy.float(minmax[3]),numpy.int(nval[1]))
print(minmax[0])

print(len(Rvec))
print(len(petitrvec))
ESCF=numpy.zeros([len(Rvec),len(petitrvec)])
RVECmesh=numpy.zeros([len(Rvec),len(petitrvec)])
petitrmesh=numpy.zeros([len(Rvec),len(petitrvec)])
nr=0
for grandr in Rvec:
	npetitr=0
	for petitr in petitrvec:
		RVECmesh[nr,npetitr] = grandr 
		petitrmesh[nr,npetitr] = petitr
		npetitr = npetitr +1
	nr = nr + 1

#        f.write(str(dOz)+" "+str(O[2]+dOz-H1[2])+" "+str(mf.energy_tot())+" "+str(mycc.energy())+'\n')
#        mc = mcscf.CASCI(mf, 8, 8)

SCFPES = numpy.load('SCFPES.npy')


# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(petitr,grandr,ESCF, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('test.png')
plt.close()
plt.plot(Rvec,SCFPES[:,5])
plt.close()
def creategrandrfig():
 for i in range(len(SCFPES[0,:])):
   plt.plot(Rvec,SCFPES[:,i])
   vec = numpy.array(SCFPES[:,i], dtype=float)
   minval=numpy.min(vec)
   print(i)
   maxval=SCFPES[-1,i]
   print(minval)
   print(maxval)
   plt.axis([Rvec[0],Rvec[-1],minval-(0.1*(maxval-minval)),maxval+0.25*(maxval-minval)])
   plt.savefig("grandr_"+str(i)+".png")
   plt.close()

def createpetitrfig():
 for i in range(len(SCFPES[:,0])):
   plt.plot(petitrvec,SCFPES[i,:])
   vec = numpy.array(SCFPES[i,:], dtype=float)
   minval=numpy.min(vec)
   maxval=SCFPES[i,-1]
   plt.axis([Rvec[0],Rvec[-1],minval-(0.1*(maxval-minval)),maxval+0.25*(maxval-minval)])
   plt.savefig("petitr_"+str(i)+".png")
   plt.close()

createpetitrfig()

plt.savefig('test2.png')
plt.close()
leveltot = plt.MaxNLocator(nbins=25).tick_values(ESCF.min(), ESCF.max())
plt.contourf(petitrmesh,RVECmesh,ESCF,levels=leveltot)
plt.savefig('test3.png')
plt.close()
