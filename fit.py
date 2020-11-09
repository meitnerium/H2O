import numpy as np
import matplotlib.pyplot as plt
def harm(x,k,x0,E0):
    return k*(x-x0)**2/2+E0
def fitharm(x,pot):
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(harm, x, pot)
    return popt,pcov


if __name__ == "__main__":
    exemple='harmH2'
    if exemple == 'harmH2':
        import configparser
        ESCF = np.load('1dimpesh2/SCFPES.npy')

        config = configparser.ConfigParser()
        config.read('1dimpesh2/config.ini')

        nval = config['GRID']['nval'].split(',')
        minmax = config['GRID']['minmax'].split(',')
        dim1vec = np.linspace(np.float(minmax[0]), np.float(minmax[1]), np.int(nval[0]))

        min=65
        max=77
        plt.plot(dim1vec[min:max+1], ESCF[min:max+1],label='Fit')

        popt, pcov = fitharm(dim1vec[min:max+1], ESCF[min:max+1])
        plt.plot(dim1vec[min:max+1], harm(dim1vec[min:max+1], *popt),label='Fit')
        print(popt)
        plt.show()