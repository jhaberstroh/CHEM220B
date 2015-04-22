import argparse
import numpy as np
import scipy.signal as scisig
import os.path
from StringIO import StringIO

def loadxyz(filename, return_comments=False, tmax=None):
    xyz = []
    comments = []
    with open(filename) as f:
        i = 0
        data = ""
        for l in f:
            if i==0:
                n_particle = int(l)
            elif i==1:
                comments.append(l)
            else:
                data = data + l
            i+=1
            if i==(n_particle+2):
                i=0
                data = StringIO(data)
                frame = np.loadtxt(data, usecols=[1,2,3])
                xyz.append(frame)
                if len(xyz) % 1000 == 0:
                    print len(xyz)
                data = ""
                if (not tmax is None) and (len(xyz) == tmax):
                    break
    xyz = np.array(xyz)
    if return_comments:
        return xyz, comments
    return xyz

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-eqfile", nargs='+', default=None)
    parser.add_argument("-savedir", default=None)
    parser.add_argument("-savename", default="St")
    parser.add_argument("-invscale", default=1., type=float, help='dividing scale factor for Xt plots')
    parser.add_argument("file", nargs='*', help='''Name of S(t) xyz files, 
        with a nomenclature where the k-value is listed as name_k_anything.xyz. String between
        the first two underscores is taken as the wavevector.''')
    parser.add_argument("-vdwsig", type=float, default=3.4, help='VDW scaling factor, default 3.4 for distance between C atoms.')
    parser.add_argument("-ylim", type=float, nargs=2, default=[-.2, .8])
    args = parser.parse_args()
    if not args.savedir is None:
        import matplotlib
        matplotlib.use("agg")
    import matplotlib.pyplot as plt
    def safesaveplot(savedir=None, name=None):
        plt.tight_layout()
        if savedir is None:
            plt.show()
        else:
            plt.savefig(savedir+'/'+name)
            plt.clf()

    if not args.eqfile is None:
        print "Running only the Equilibrium analysis"
        S_kt_collection = []
        for eqfile in args.eqfile:
            xyz, comment = loadxyz(eqfile, return_comments=True)
            params = np.zeros((len(comment), 4))
            for j, l in enumerate(comment):
                l = l.strip()
                l = [w.split()[0] for w in l.split('=')]
                l = [float(n) for n in l[1:]]
                params[j, :] = np.array(l)
            equil = (params[:, 0] >= 0)
            # Factor included for system size scaling for C atom VDW radius
            L = np.power(250, 1./3.) * args.vdwsig
            k = np.array((12. * np.pi / L, 
                          24. * np.pi / L, 
                          36. * np.pi / L))
            x_ti = xyz[equil, :, 0]
            phase_kti = k[:, np.newaxis, np.newaxis] * x_ti[np.newaxis, :] 
            p_plus_kt = np.sum( np.cos( phase_kti ) + 1j * np.sin( phase_kti ), axis=2)
            p_minu_kt = np.sum( np.cos(-phase_kti ) + 1j * np.sin(-phase_kti ), axis=2)
            # Skt1  = scisig.fftconvolve(p_plus_kt[0, :], p_minu_kt[0, ::-1])
            # Skt2  = scisig.fftconvolve(p_plus_kt[1, :], p_minu_kt[1, ::-1])
            # Skt3  = scisig.fftconvolve(p_plus_kt[2, :], p_minu_kt[2, ::-1])
            Nt   = x_ti.shape[0]
            S_kt = np.zeros( (len(k), (2*Nt)-1), dtype=np.complex128 )
            S_kt[0, :] = scisig.fftconvolve(p_plus_kt[0, :], p_minu_kt[0, ::-1])
            S_kt[1, :] = scisig.fftconvolve(p_plus_kt[1, :], p_minu_kt[1, ::-1])
            S_kt[2, :] = scisig.fftconvolve(p_plus_kt[2, :], p_minu_kt[2, ::-1])
            S_kt = S_kt[:, Nt:] / np.linspace(Nt, 1, Nt-1) / x_ti.shape[1]
            t_t =  params[equil, 0][:-1]
            sub = (t_t < 1.0)
            t_t  = t_t[sub]
            S_kt = S_kt[:,sub]
            S_kt_collection.append(S_kt / 2.)
        S_kt_collection = np.array(S_kt_collection)
        S_kt_mean = np.mean(S_kt_collection, axis=0)
        S_kt_rerr  = np.std(S_kt_collection.real, ddof=1, axis=0) \
                    / np.sqrt(S_kt_collection.shape[0] - 1)
        S_kt_ierr  = np.std(S_kt_collection.imag, ddof=1, axis=0) \
                    / np.sqrt(S_kt_collection.shape[0] - 1)
        print S_kt_ierr
        for i in range(S_kt_mean.shape[0]):
            plt.plot(        t_t, S_kt_mean[i, :].real)
            plt.fill_between(t_t, S_kt_mean[i, :].real - S_kt_rerr[i, :],
                                  S_kt_mean[i, :].real + S_kt_rerr[i, :],
                             alpha = .2, facecolor='blue')
        for i in range(S_kt_mean.shape[0]):
            plt.plot(        t_t, S_kt_mean[i, :].imag)
            plt.fill_between(t_t, S_kt_mean[i, :].imag - S_kt_ierr[i, :],
                                  S_kt_mean[i, :].imag + S_kt_ierr[i, :],
                             alpha = .2, facecolor='red')
        plt.legend(['k=12, REAL',
                    'k=24, REAL',
                    'k=36, REAL',
                    'k=12, IMAG',
                    'k=24, IMAG',
                    'k=36, IMAG'])
        plt.ylim(args.ylim)
        plt.title("S(t) estimate, with error bars")
        safesaveplot(args.savedir, 'St_all.png')

    else: 
        print "Running only the dynamical analysis"
        Xt_coll = []
        Xt_k    = []
        for i, file in enumerate(args.file):
            print file
            xyz, comment = loadxyz(file, return_comments=True)
            file_k = os.path.basename(file).split("_")[1]
            Xt_k.append(file_k)

            params = np.zeros((len(comment), 4))
            for j, l in enumerate(comment):
                l = l.strip()
                l = [w.split()[0] for w in l.split('=')]
                l = [float(n) for n in l[1:]]
                params[j, :] = np.array(l)
            relax = (params[:,3] == 0.0)
            L = params[0,1]
            k = params[0,2]
            print k, L
            T = params[relax, 0]
            x_ti  = xyz[relax, :, 0] / args.vdwsig
            x_tij = x_ti[:, np.newaxis,  :] - x_ti[0, :, np.newaxis] 
            Xt_fourier = np.sum( np.cos( x_ti  * k ), axis = 1) / np.sqrt(x_ti.shape[1])
            Xt_coll.append(Xt_fourier)
            # plt.plot(params[relax, 0], fourier)
            # safesaveplot(args.savedir, '{}_{}.png'.format(args.savename, i))
        k_vals  = set(Xt_k)
        Xt_k    = np.array(Xt_k)
        Xt_coll = np.array(Xt_coll).T
        print "K vals:", k_vals
        legend = []
        for k in k_vals:
            sub = (Xt_k == k)
            print Xt_coll.shape
            Xt_coll_k = Xt_coll[:, sub] / args.invscale
            Xt_mean = np.mean(Xt_coll_k, axis=1)
            Xt_yerr = np.std(Xt_coll_k, axis=1, ddof=1) / np.sqrt(Xt_coll_k.shape[1] - 1)
            plt.plot(params[relax, 0], Xt_mean)
            plt.fill_between(params[relax, 0], Xt_mean - Xt_yerr, Xt_mean + Xt_yerr,
                    alpha = .2, facecolor='blue')
            legend.append("k={}".format(k))
        plt.ylim(args.ylim)
        plt.legend(legend)
        plt.title("Mean X(t) for Bf={}, with error estimates".format(args.invscale))
        safesaveplot(args.savedir, '{}.png'.format(args.savename))

        # plt.plot(params[relax, 0], Xt_coll)
        # plt.title("Collection of X(t) values")
        # safesaveplot(args.savedir, '{}_XtColl.png'.format(args.savename))
    



if __name__=="__main__":
    main()
