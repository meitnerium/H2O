
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





nval=config['GRID']['nval'].split(',')
minmax=config['GRID']['minmax'].split(',')
dim1vec=numpy.linspace(numpy.float(minmax[0]),numpy.float(minmax[1]),numpy.int(nval[0]))
ESCF=numpy.load('SCFPES.npy')
MO1=numpy.load('MO1.npy')
MO2=numpy.load('MO2.npy')

plt.plot(dim1vec,MO1)
plt.plot(dim1vec,MO2)




#plt.plot(dim1vec,ESCF)
minval=numpy.min(ESCF)
maxval=ESCF[-1]
#plt.axis([dim1vec[0],dim1vec[-1],minval-(0.1*(maxval-minval)),maxval+0.25*(maxval-minval)])


plt.savefig('test.png')
plt.show()
plt.close()
