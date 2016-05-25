import numpy as N
import matplotlib.pyplot as plt
import scipy.interpolate as I



def xilin_from_pklin(r, pklin=None, smw=3.0, L=20000, upper=4*N.pi, lower=10**-3):
    k=pklin[:,0]
    pk=pklin[:,1]
    pklin=I.interp1d(k, pk)
    if (upper > k.max()):
        upper=k.max()
    deltak=upper/L
    kk=N.arange(lower, deltak*L, deltak)
    sk=r.shape
    rkk=N.outer(r,kk)
    jrk=N.sin(rkk)/rkk
    jrk[N.where(rkk==0)]=1.0
    w=N.exp(-smw**2*kk*kk/2)
    ans=N.sum((pklin(kk)*w**2*kk**2*jrk*deltak)/(2*N.pi**2), 1)
    ans=N.reshape(ans, sk)
    return ans



def plot_gaussian(pkfile="/gpfs/data/nmccull/rvb/pkinit.txt", smoothing=3.0):
    r=N.arange(1.0, 200, 1.0)	
    pklin=N.loadtxt(pkfile)
    xilin=xilin_from_pklin(r, pklin=pklin, smw=smoothing)
    delta=N.linspace(-5, 5, 50)
    sigma2=xilin_from_pklin(N.array([0.0]), pklin=pklin, smw=smoothing)
    #this is xi(r, delta) (but not multiplied by delta)
    xird=N.outer(delta*N.exp(-delta**2/(2*sigma2)), xilin/N.sqrt(2*N.pi)/sigma2**(3./2.))
    rmat=N.tile(r, (delta.size, 1))
    dmat=N.transpose(N.tile(delta, (r.size, 1)))
    plt.pcolor(r, delta, N.arcsinh(xird*rmat**2*dmat), cmap=plt.get_cmap('gist_earth'))
    plt.colorbar()
    plt.ylim([-5, 5])
    plt.xlim([1.0, 130])
    plt.xlabel('$r$ [Mpc/$h$]')
    plt.ylabel(r'$\delta$')
    plt.savefig('/Users/nuala/Documents/Research/Figures/slicedcorr/model_gauss.png', bbox_inches='tight')
    plt.show()

def plot_lognormal(pkfile="/gpfs/data/nmccull/rvb/pkinit.txt", smoothing=3.0):
    r=N.arange(1.0, 200, 1.0)	
    pklin=N.loadtxt(pkfile)
    xilin=xilin_from_pklin(r, pklin=pklin, smw=smoothing)
    logd=N.linspace(-5, 5, 50)
    sigma2=xilin_from_pklin(N.array([0.0]), pklin=pklin, smw=smoothing)
    mu=-sigma2/2.0
    print sigma2, mu
    
    Amat=N.transpose(N.tile(logd, (r.size, 1)))
    ximat=N.tile(xilin, (logd.size, 1))
    rmat=N.tile(r, (logd.size, 1))
    t1=N.exp((sigma2**2+mu*(sigma2-ximat)+(Amat-ximat)*ximat)**2/(2*sigma2*(sigma2**2-ximat**2)))
    t2=N.exp((mu*(sigma2-ximat)+Amat*ximat)**2/(2*sigma2*(sigma2**2-ximat**2)))
    t3=N.exp(-(Amat**2*sigma2+2*mu**2*(sigma2-ximat)+2*Amat*mu*(-sigma2+ximat))/(2*(sigma2**2-ximat**2)))

    #this is xi(r, delta) (but not multiplied by delta)
    xird=t3*(t1-t2)/N.sqrt(2*N.pi*sigma2)

    plt.pcolor(r, logd, N.arcsinh(xird*rmat**2*(N.exp(Amat)-1)), cmap=plt.get_cmap('gist_earth'))
    plt.colorbar()
    plt.ylim([-5, 5])
    plt.xlim([1.0, 130])
    plt.xlabel('$r$ [Mpc/$h$]')
    plt.ylabel(r'$\ln(1+\delta)$')
    plt.savefig('/Users/nuala/Documents/Research/Figures/slicedcorr/model_lognormal.png', bbox_inches='tight')
    plt.show()
