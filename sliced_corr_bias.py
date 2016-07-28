import numpy as N
import scipy.interpolate as I
import scipy.integrate as INT
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

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


# computes the halo density for a given matter density (array) given the parameters (b, alpha, epsilon, rhoe), with or without exponential cutoff
# (if include_exp=0, epsilon and rhoe are ignored)
def delta_h(dm, b, alpha, epsilon, rhoe):
    rho_m=dm+1.0
    am=N.log(rho_m)
    
    ah=N.log(b)+alpha*am
    if (rhoe!=0):
        ae=N.log(rhoe)
        term=-N.exp(-epsilon*(am-ae))
        ah=ah+term
    return N.exp(ah)-1.0
 
# returns the derivative of delta_h with respect to delta_m in a power law bias model (with exponential cutoff if include_exp=1)   
def dhdm(dm, b, alpha, epsilon, rhoe):
    rho_m=dm+1
    if (rhoe!=0):
        return alpha*b*rho_m**(alpha-1.0)*N.exp(-(rho_m/rhoe)**-epsilon)+b*epsilon*rho_m**(alpha-epsilon-1.0)*N.exp(-(rho_m/rhoe)**-epsilon)*rhoe**epsilon
    else:
        return alpha*b*rho_m**(alpha-1.0)
 
# returns the derivative of log 1+delta_h with respect to log 1+delta_m in a power law bias model (with exponential cutoff if include_exp=1)
def dahdam(dm, b, alpha, epsilon, rhoe):
    rho_m=dm+1
    am=N.log(rho_m)
    if (rhoe!=0):
        ae=N.log(rhoe)
        return alpha+epsilon*N.exp(-epsilon*(am-ae))
    else:
        return alpha*N.ones_like(am)
    
# does the inverse transform of delta_h = g(delta_m)
# so this returns delta_m = g^-1(delta_h)
# it does this using interpolation
def delta_m(dh, bias_func, b, args):
    dm2=N.arange(-1.0, 5*dh.max(), 0.01)
    
    dh2=bias_func(dm2, b, *args)
    i=1
    while (dh2.max()<dh.max()):
        dm2=N.arange(-1.0, (5+i)*dh.max(), 0.01)
        dh2=bias_func(dm2, b, *args)
        i=i+1
    
    
    dmi=I.interp1d(dh2, dm2)
    return dmi(dh)
    
    


def bigauss(x1, x2, xilin, sigma):
    rho=xilin/sigma**2
    return 1./(2*N.pi*sigma**2*N.sqrt(1-rho**2))*N.exp(-(x1**2+x2**2-2*rho*x1*x2)/(2*sigma**2*(1-rho**2)))
    

def integrand(am1, *args):
    am2=args[0]
    xilin=args[1]
    sigma=args[2]
    b=args[3]
    alpha=args[4]
    epsilon=args[5]
    rhoe=args[6]
    mu=args[7]
    
    dm1=N.exp(am1)-1.0
    dh1=delta_h(dm1, b, alpha, epsilon, rhoe)
    return bigauss(am1-mu, am2-mu, xilin, sigma)*dh1
    
    
def gauss(x1, mu, sigma):
    return 1./(sigma*N.sqrt(2*N.pi))*N.exp(-(x1-mu)**2/(2*sigma**2))


# This is the main function in the program
# dh0 is an array of the values of delta_h at which you want to compute the sliced correlation function
# xilin is an array containing the (matter) correlation function (<log rho1 log rho2>) at the values of r that you want the sliced correlation function
# sigma is the standard deviation of the log rho field (matter density)
# alpha, epsilon, and rhoe are parameters of the bias model (given in the function delta_h)
# if rhoe=0, the exponential part is not included
# bins='log' returns xi(r, delta) using logarithmic binning in delta
# otherwise it returns it in linear delta bins
# returns the halo density (or log density if bins='log') and the sliced correlation function
def compute_xird(dh0, xilin, sigma, bias_func, args, bins='log'):
    mu=-sigma**2/2.0
    
    amin=-5*sigma
    amax=5*sigma
    
    #figure out b by computing the mean of delta_h
    a0=N.linspace(amin, amax, 200)
    da=N.diff(a0)[0]
    fg=gauss(a0, mu, sigma)
    
    rhom=N.exp(a0)
    dh=bias_func(rhom-1, 1.0, *args)
    mean=N.sum(fg*(1+dh)*da)
    b=1.0/mean
    
    xird=N.zeros((dh0.size, xilin.size), dtype=N.float)
    dm0=delta_m(dh0, bias_func, b, args)
    aa0=N.log(1+dm0)
    
    ah0=N.log(1+dh0)
    return ah0, aa0
    print b
    
    if (bins=='log'):
        dmt=dm0[:-1]+N.diff(dm0)/2.0
        dahdam=InterpolatedUnivariateSpline(dmt, N.diff(ah0)/N.diff(aa0), k=1)
        factor=1.0/dahdam(dm0)
    else:
        dmt=dm0[:-1]+N.diff(dm0)/2.0
        dhdm=InterpolatedUnivariateSpline(dmt, N.diff(dh0)/N.diff(dm0), k=1)
        factor=(1.0)/(1.0+dm0)*1.0/dhdm(dm0)
    
    return dm0, dahdam(dm0)
    for i in N.arange(dh0.size):
        for j in N.arange(xilin.size):
            xird[i, j]=INT.quad(integrand, amin, amax, args=(aa0[i], xilin[j], sigma, b, alpha, epsilon, rhoe, mu))[0]
            xird[i,j]=xird[i,j]*dh0[i]*factor[i]
    if (bins=='log'):
        return N.log(1+dh0), xird
    else:
        return dh0, xird
    
    
    
# This function shows how to call the above routine and plot it
def plot_xird(alpha=1.0, epsilon=1.5, rhoe=0.4, pkfile="/gpfs/data/nmccull/rvb/pkinit.txt", figfile="", smoothing=3.0):

    rr=N.arange(1.0, 200.0, 1.0)
    pklin=N.loadtxt(pkfile)
    xilin=xilin_from_pklin(rr, pklin=pklin, smw=smoothing)
    sigma=N.sqrt(xilin_from_pklin(N.array([0.0]), pklin=pklin, smw=smoothing))
    logd=N.linspace(-5, 5, 50)
    
    a0, xird=compute_xird(N.exp(logd)-1.0, xilin, sigma, alpha, epsilon,rhoe, bins='log')
    
    plt.pcolor(rr, a0, N.arcsinh(rr**2*xird), cmap='gist_earth')
    plt.xlim([5, 130])
    plt.ylim([-5, 5])
    plt.colorbar()
    plt.xlabel('r [Mpc/$h$]')
    plt.ylabel(r'$\ln(1+\delta_h)$')
    if (figfile!=""):
        plt.savefig(figfile, bbox_inches='tight')
    plt.show()
    
    
    