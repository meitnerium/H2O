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

def numerovwindows(pot,Emax,Emin,dx,tol,nodetofound):
    psi=np.zeros(len(pot))
    n=1
    while Emax-Emin > tol:
        E=Emin+(Emax-Emin)/2
        jmin=3
        print(n,Emin,E,Emax)
        psi=np.zeros(len(pot))
        psi[0] = 0
        psi[1] = 1E-4

        negderfound=False
        for i in range(len(pot)-2):
            j=i+2
            psi[j]=numerovplus(psi[j-2],psi[j-1],pot[j-2],pot[j-1],pot[j],E,dx)
            der=(psi[j]-psi[j-1])/dx
            if der < 0:
                jmin=j
                print(jmin)
                negderfound=True
                break
        psi[0:jmin] = psi[0:jmin] / psi[jmin]
        der=(psi[j]-psi[j-1])/dx
        psi[-1]=0
        if nodetofound%2 == 0:
            psi[-2]=1E-4
        else:
            psi[-2]=-1E-4
        psi[-2]=1E-4
        node=0
        maxfound=False
        if negderfound:
          for j in range(len(pot)-3,jmin-1,-1):
            psi[j]=numerovplus(psi[j+2],psi[j+1],pot[j+2],pot[j+1],pot[j],E,dx)
            if psi[j]*psi[j+1] < 0:
                print("node found")
                node=node+1
                if node>nodetofound:
                    maxfound=True
                    break
          psi[jmin:] = psi[jmin:] / psi[jmin]
          print('der1')
          print(der)
          print('der2')
          der2=(psi[jmin+1]-psi[jmin])/dx
          print((psi[jmin+1]-psi[jmin])/dx)
          if node < nodetofound:
            print('Energy too low')
            Emin=E
          elif der*der2>0 and maxfound==False:
            print('Energy too low')
            Emin=E
          else:
            Emax=E

          psi[jmin:]=psi[jmin:]/psi[jmin]*psi[jmin-1]
          print(psi[jmin-1])
          print(psi[jmin])
          print(psi[jmin+1])
        else:
          Emin=E
        n=n+1

        #plt.plot(psi)
        #plt.show()
        #plt.savefig('test.png')
        #plt.close()
        #exit(1)
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

        psi[jmin:]=psi[jmin:]/psi[jmin]*psi[jmin-1]
    print('nodetobefound = '+str(nodetofound))
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
    exemple='H2plus'
    if exemple == 'harm':
        numpoints=1024
        xmin=-5
        xmax=5
        x=np.linspace(xmin,xmax,numpoints)
        pot=0.5*x**2
        E,psi = numerovwindows(pot, 10, 0, x[1]-x[0],1E-10,0)
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

        E,psi = numerovwindows(ESCF, 0, -3, dim1vec[1] - dim1vec[0], 1E-10,0)
        plt.plot(dim1vec,psi)
        np.save('1dimpesh2/h2fund',psi)
        plt.savefig('test.png')
        #plt.show()
    elif exemple == 'H2plus':
        # H2+
        maxv=12
        import configparser
        ESCF = np.load('1dimpesh2plus/SCFPES.npy')

        config = configparser.ConfigParser()
        config.read('1dimpesh2plus/config.ini')

        nval = config['GRID']['nval'].split(',')
        minmax = config['GRID']['minmax'].split(',')
        print('min='+str(minmax[0]))
        print('max=' + str(minmax[1]))

        dim1vec = np.linspace(np.float(minmax[0]), np.float(minmax[1]), np.int(nval[0]))
        allpsi=[]
        v=0
        while v <= maxv: 
            E,psi = numerovwindows(ESCF[0:257], 3, -3, dim1vec[1] - dim1vec[0], 1E-10,v)
            #plt.plot(dim1vec[0:257],psi)
            plt.plot(psi)
            plt.show()
            v=v+1
        plt.savefig('test.png')
        #plt.show()
    elif exemple == 'test':
        import configparser
        ESCF = np.load('1dimpesh2plus/SCFPES.npy')

        config = configparser.ConfigParser()
        config.read('1dimpesh2plus/config.ini')

        nval = config['GRID']['nval'].split(',')
        minmax = config['GRID']['minmax'].split(',')
        print(nval)
        dim1vec = np.linspace(np.float(minmax[0]), np.float(minmax[1]), np.int(nval[0]))
        allpsi=[]
        v=10
        Emin=-8.931149137020111e-3
        Emax=0.1
        E,psi = numerovwindows(ESCF[0:257], Emax, Emin, dim1vec[1] - dim1vec[0], 1E-15,v)
        plt.plot(dim1vec,psi)
        plt.show()
        plt.savefig('test.png')
        #plt.show()

