import numpy as np
import matplotlib.pyplot as plt
def numerovplus(psimoin2,psimoin1,potmoins2,potmoins1,pot,E,dx):
    G1=2*(potmoins2-E)
    G2=2*(potmoins1-E)
    G3=2*(pot-E)
    num=psimoin1*2-psimoin2+5*G2*psimoin1*dx**2/6+G1*psimoin2*dx**2/12
    denom=1-G3*dx**2/12
    psi=num/denom
    return psi

def numerovwindows(pot,Emax,Emin,dx,tol):
    psi=np.zeros(len(pot))
    while Emax-Emin > tol:
        E=Emin+(Emax-Emin)/2
        psi=np.zeros(len(pot))
        psi[0] = 0
        psi[1] = 1E-4

        for i in range(len(pot)-2):
            j=i+2
            psi[j]=numerovplus(psi[j-2],psi[j-1],pot[j-2],pot[j-1],pot[j],E,dx)
            der=(psi[j]-psi[j-1])/dx
            if der < 0:
                jmin=j
                break

        psi[-1]=0
        psi[-2]=1E-4
        for j in range(len(pot)-3,i-1,-1):
            psi[j]=numerovplus(psi[j+2],psi[j+1],pot[j+2],pot[j+1],pot[j],E,dx)
            if psi[j]*psi[j+1] < 0:
                print("node found")
                Emax=E
                break

    plt.plot(psi)
    plt.savefig('test.png')
    plt.close()
    exit(1)
        #test,i,psi = numerov(pot,E,dx)
        #print(test)
        #if test == 1:
        #    print("E too big")
        #    Emax=E
        #if test == 2:
        #    numerovinv(psi,pot,E,dx,i)
        #    break
        #else:
        #    print("E too low")
        #    Emin=E

    print('E final = '+str(E))
    return E,psi

def numerovinv(psi,pot,E,dx,i):
    print(i)    
    psi[-1]=0
    psi[-2]=1E-4
    
    for j in range(len(pot)-3,i-1,-1):
        psi[j]=numerovplus(psi[j+2],psi[j+1],pot[j+2],pot[j+1],pot[j],E,dx)
    
def numerov(pot,E,dx):
    psi=np.zeros(len(pot))

    print("E="+str(E))
    psi[0] = 0
    psi[1] = 1E-4
    for i in range(len(pot)-2):
        j=i+2
        psi[j]=numerovplus(psi[j-2],psi[j-1],pot[j-2],pot[j-1],pot[j],E,dx)
        if psi[j] < 0:
            print("neg function")
            return 1,i,psi
        der=(psi[j]-psi[j-1])/dx
        if der < 0:
            print("neg derivative")
            return 2,i,psi
    #plt.plot(psi)
    #plt.show()
    return -1,i,psi

if __name__ == "__main__":
    exemple='H2'
    if exemple == 'harm':
        numpoints=1024
        xmin=-5
        xmax=5
        x=np.linspace(xmin,xmax,numpoints)
        pot=0.5*x**2
        E,psi = numerovwindows(pot, 10, 0, x[1]-x[0],1E-10)
        plt.plot(x,psi)
        plt.savefig('test.png')
        #plt.show()

    elif exemple == 'H2':
        # H2
        import configparser
        ESCF = np.load('1dimpesh2/SCFPES.npy')

        config = configparser.ConfigParser()
        config.read('1dimpesh2/config.ini')

        nval = config['GRID']['nval'].split(',')
        minmax = config['GRID']['minmax'].split(',')
        dim1vec = np.linspace(np.float(minmax[0]), np.float(minmax[1]), np.int(nval[0]))

        E,psi = numerovwindows(ESCF, 0, -3, dim1vec[1] - dim1vec[0], 1E-10)
        plt.plot(dim1vec,psi)
        np.save(psi,'1dimpesh2/h2fund')
        plt.savefig('test.png')
        #plt.show()
    elif exemple == 'H2plus':

        # H2+
        import configparser
        ESCF = np.load('1dimpesh2plus/SCFPES.npy')

        config = configparser.ConfigParser()
        config.read('1dimpesh2plus/config.ini')

        nval = config['GRID']['nval'].split(',')
        minmax = config['GRID']['minmax'].split(',')
        dim1vec = np.linspace(np.float(minmax[0]), np.float(minmax[1]), np.int(nval[0]))

        E,psi = numerovwindows(ESCF, 0, -3, dim1vec[1] - dim1vec[0], 1E-10)
        plt.plot(dim1vec,psi)
        plt.savefig('test.png')
        #plt.show()
