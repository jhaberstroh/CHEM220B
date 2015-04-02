import numpy as np
import argparse

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()

parser = argparse.ArgumentParser()
parser.add_argument("-NK", type=int, default=1000, help="Number of K points to use")
parser.add_argument("-save", help="Number of K points to use")
args = parser.parse_args()

if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

def percus_yevick():
    show_basis = False
    
    K = np.linspace(np.sqrt(.018), np.sqrt(100), args.NK)
    K = np.square(K)
    C1 = - np.cos(K)/np.power(K, 2) + np.sin(K)/np.power(K, 3)
    C2 = - np.cos(K)/np.power(K, 2) + 2. * np.sin(K)/np.power(K, 3) +\
            2. * (np.cos(K) - 1.)/np.power(K, 4)
    C3 = - np.cos(K)/np.power(K, 2) + 4. * np.sin(K)/np.power(K, 3) +\
            12. * np.cos(K)/np.power(K, 4) - 24. * np.sin(K)/np.power(K, 5) -\
            24. * (np.cos(K) - 1)/np.power(K, 6)
    K[0] = 0
    C1[0] = 1./3.
    C2[0] = 1./4.
    C3[0] = 1./6.
    
    if show_basis:
        plt.plot(K, C1)
        plt.plot(K, C2)
        plt.plot(K, C3)
        plt.hlines([1./6., 1./4., 1./3.], 0, 10)
        safesaveplot(args.save, "basisfn.png")
    
    densities = [.02, .1, .2, .3, .4, .5, .6, .7, .8, .9]
    for density in densities:
        packing = (np.pi / 6.) * density
        lam1 = -(1. + 2. * packing)**2/(1 - packing)**4
        lam2 = 6*packing * (1. + packing/2.)**2 / (1 - packing)**4
        lam3 = packing / 2. * lam1
        Ctot_k = (lam1 * C1 + 
                  lam2 * C2 + 
                  lam3 * C3) * np.pi * 4.
        
        I_k = np.square(Ctot_k) / (1 - density * Ctot_k)
        
        Nr = 1000
        r = np.linspace(.01, 10, Nr)
        C_r = lam1 + lam2 * r + lam3 * np.power(r, 3)
        C_r[r >= 1.] = 0
        I_r = np.zeros(Nr)
        dK = np.diff(K)
        for i, r_i in enumerate(r):
            integrand = I_k * np.sin(K * r_i) * K
            value = np.sum((integrand[1:] + integrand[:-1]) / 2. * dK)
            I_r[i] = value * density / (2 * np.pi**2 * r_i)
        
        sub = r >= 0
        plt.plot(r[sub], I_r[sub]+C_r[sub], c=(max(density*2-1, 0), 0, max(1-density*2., 0)))
    plt.ylim([-1.3,4.5])
    plt.xlim([0,5])
    fig = plt.gcf()
    fig.set_size_inches(6.0, 6.5)
    plt.legend(densities)
    plt.title("PY h(r)".format(density))
    safesaveplot(args.save, "percus_yevick_h.png")

    densities = [.1, .3, .7, .9]
    for density in densities:
        packing = (np.pi / 6.) * density
        lam1 = -(1. + 2. * packing)**2/(1 - packing)**4
        lam2 = 6*packing * (1. + packing/2.)**2 / (1 - packing)**4
        lam3 = packing / 2. * lam1
        Ctot_k = (lam1 * C1 + 
                  lam2 * C2 + 
                  lam3 * C3) * np.pi * 4.
        
        I_k = np.square(Ctot_k) / (1 - density * Ctot_k)
        
        Nr = 1000
        r = np.linspace(.01, 10, Nr)
        C_r = lam1 + lam2 * r + lam3 * np.power(r, 3)
        C_r[r >= 1.] = 0
        I_r = np.zeros(Nr)
        dK = np.diff(K)
        for i, r_i in enumerate(r):
            integrand = I_k * np.sin(K * r_i) * K
            value = np.sum((integrand[1:] + integrand[:-1]) / 2. * dK)
            I_r[i] = value * density / (2 * np.pi**2 * r_i)
        
        sub = r >= 1
        plt.plot(r[sub], 1+I_r[sub]+C_r[sub])
        plt.ylim([0,4.5])
        plt.xlim([0,5])
        fig = plt.gcf()
        fig.set_size_inches(3.0, 2.5)
        d = int(round(density * 100))
        plt.title("PY g(r), density = {}".format(density))
        safesaveplot(args.save, "percus_yevick_g{}.png".format(d))


def harmonic_md_plot(dt, Tf = 6.5):
    Nt = int(np.ceil(Tf / dt))
    x = np.zeros(Nt)
    v = np.zeros(Nt)
    T = np.array(range(Nt)) * dt
    x[0] = 1.
    v[0] = 0.
    for i, t in enumerate(T):
        if i != 0:
            v[i] = v[i-1] - .5 * x[i-1] * dt
            x[i] = x[i-1] +      v[ i ] * dt
            v[i] = v[ i ] - .5 * x[ i ] * dt
    print x.shape, T.shape
    plt.plot(T, x)
    theta = np.arctan2( dt * np.sqrt(4. - dt**2), 2. - dt**2)
    plt.plot(T, np.cos(T*theta / dt), 'o')
    plt.legend(["Simulated", "Predicted"])
    plt.xlim([0,Tf])
    safesaveplot(args.save, "md-traj-{}.png".format(dt))

def harmonic_md():
    fig = plt.gcf()
    fig.set_size_inches(3.0, 3.0)
    harmonic_md_plot(.1)
    harmonic_md_plot(.4)
    harmonic_md_plot(.8)
    harmonic_md_plot(1.5)
    harmonic_md_plot(2)
    harmonic_md_plot(2.1, Tf=20.)

    fig = plt.gcf()
    fig.set_size_inches(6.0, 6.5)
    dt = .1
    Nt = 100
    x = np.zeros(Nt)
    v = np.zeros(Nt)
    T = np.array(range(Nt)) * dt
    x[0] = 1.
    v[0] = 0.
    for i, t in enumerate(T):
        if i != 0:
            v[i] = v[i-1] - .5 * x[i-1] * dt
            x[i] = x[i-1] +      v[ i ] * dt
            v[i] = v[ i ] - .5 * x[ i ] * dt
    plt.plot(T, .5 * np.square(x))
    plt.plot(T, .5 * np.square(v))
    plt.plot(T, .5 * (np.square(x) + np.square(v)))
    safesaveplot(args.save, "md-energy.png")

    plt.plot(v, x)
    safesaveplot(args.save, "md-phase.png")


    dt = .1
    Nt = 10000
    x = np.zeros(Nt)
    v = np.zeros(Nt)
    T = np.array(range(Nt)) * dt
    x[0] = 1.
    v[0] = 0.
    for i, t in enumerate(T):
        if i != 0:
            v[i] = v[i-1] - .5 * x[i-1] * dt
            x[i] = x[i-1] +      v[ i ] * dt
            v[i] = v[ i ] - .5 * x[ i ] * dt
    hist, bins = np.histogram(x, density=True, bins=500)
    plt.plot( (bins[1:] + bins[:-1])/2. , hist)
    plt.xlim([-1.1, 1.1])
    safesaveplot(args.save, "md-probability.png")

percus_yevick()
fig = plt.gcf()
fig.set_size_inches(6.0, 6.5)
harmonic_md()
