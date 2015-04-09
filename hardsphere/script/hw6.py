import numpy as np
import numpy.fft as FFT
import scipy.signal
import argparse
from StringIO import StringIO

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument("-velfile", nargs='+', required=True, help="Path to velocity data")
parser.add_argument("-xyzfile", nargs='+', required=True, help="Path to xyz data")
parser.add_argument("-enerfile", required=True, help="Path to energy data")
parser.add_argument("-save", help="Path to save images out")
parser.add_argument("-savename", default="dx", help="Base name for output")
parser.add_argument("-density", required=True, type=float, help="Density of simulation, in N/V")
args = parser.parse_args()
if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt


def loadxyz(filename):
    xyz = []
    with open(filename) as f:
        i = 0
        data = ""
        for l in f:
            if i==0:
                n_particle = int(l)
            elif i==1:
                comment = l
            else:
                data = data + l
            i+=1
            if i==(n_particle+2):
                i=0
                data = StringIO(data)
                frame = np.loadtxt(data, usecols=[1,2,3])
                xyz.append(frame)
                data = ""
    xyz = np.array(xyz)
    return xyz
        
def pos_plot():
    dt_pos = .1
    D_xx_set = []
    for i, xyzfile in enumerate(args.xyzfile):
        with open(xyzfile) as f:
            N = int(f.readline())
        xyz = loadxyz(xyzfile)
        L = np.power(N / args.density, 1./3.)
        xyz = loadxyz(xyzfile)
        shift = L/2. - xyz[0,:,:] 
        xyz = (xyz + shift + L) % L
        xyz -= L/2.
        xyz = xyz.reshape(xyz.shape[0], xyz.shape[1] * 3)
        xyz_diff = np.diff(xyz, axis=0)
        xyz_diff[xyz_diff >  L/2.] -= L
        xyz_diff[xyz_diff < -L/2.] += L
        xyz = np.cumsum(xyz_diff, axis=0)
        t = np.array(range(xyz.shape[0])) * dt_pos
        plt.plot(t, xyz)
        plt.title("Ensemble for rho={}".format(args.density))
        plt.xlabel("Time, tau")
        safesaveplot(args.save, "{}_{}_x_ensemble.png".format(args.savename, i))

        rmsd = np.std(xyz, axis=1)
        t = np.array(range(1, len(rmsd)+1)) * dt_pos
        plt.plot(t, rmsd)
        plt.title("Dynamical RMSD, rho={}".format(args.density))
        plt.xlabel("Time, tau")
        safesaveplot(args.save, "{}_{}_x_rmsd.png".format(args.savename, i))

        diff = np.square(rmsd) / t / 2.
        plt.plot(t, diff)
        safesaveplot(args.save, "{}_{}_x_diffusion.png".format(args.savename, i))

        D_xx = np.mean(diff[-10:])
        D_xx_set.append(D_xx)
    if len(D_xx_set) > 1:
        D_xx_set = np.array(D_xx_set)
        D_xx_est = np.mean(D_xx_set)
        D_xx_err = np.std(D_xx_set, ddof=1) / np.sqrt(len(D_xx_set))
        print "D_xx = {} +/- {}".format(D_xx_est, D_xx_err)
    else:
        print "D_xx = {}".format(D_xx_set[0])

def vel_plot():
    dt_vel = .01
    D_vv_set = []
    for fnum, velfile in enumerate(args.velfile):
        # Velocities are indexed with [T, 3i + dir],
        #  where i indexes the particle
        vel = np.loadtxt(velfile, skiprows=1)
        vel = vel.reshape((vel.shape[0], vel.shape[1]/3, 3))
        t = np.array(range(vel.shape[0])) * dt_vel
        plt.plot(t, vel[:,0,0])
        safesaveplot(args.save, "{}_{}_v_particle.png".format(args.savename, fnum))

        Ct = np.zeros((vel.shape[0] - 1))
        t = np.array(range(len(Ct))) * dt_vel
        for i in xrange(vel.shape[1]):
            for j in xrange(vel.shape[2]):
                sig = vel[:,i,j]
                Ct_conv = scipy.signal.fftconvolve(sig, sig[::-1])
                T = len(Ct_conv) / 2
                Ct_conv = Ct_conv[T+1:]
                norms = T - np.array(range( T ))
                Ct += Ct_conv / norms
        Ct /= (vel.shape[1] * vel.shape[2])
        plt.plot(t, Ct)
        plt.title("Velocity Autocovariance, rho={}".format(args.density))
        plt.xlabel("Time, tau")
        safesaveplot(args.save, "{}_{}_v_autocov.png".format(args.savename, fnum))
        
        plt.hist(vel.flatten(), log=True, bins=100, normed=True)
        kT = 1.5
        sig = kT / 2.
        prob_curve_x = np.linspace(-5, 5, 100)
        prob_curve_y = 1. / (np.sqrt(2. * np.pi) * sig) * \
                    np.exp(-np.square(prob_curve_x)/(2. * sig))
        plt.plot(prob_curve_x, prob_curve_y, 'r--')
        plt.title("Velocity distribution, one component, rho={}".format(args.density))
        plt.xlabel("Velocity, sigma/tau")
        plt.ylabel("Probability")
        safesaveplot(args.save, "{}_{}_v_hist.png".format(args.savename, fnum))

        tf = 350
        D_vv = np.sum((Ct[1:tf] + Ct[:tf-1]) / 2.0) * dt_vel
        D_vv_set.append(D_vv)
    if len(D_vv_set) > 1:
        D_vv_set = np.array(D_vv_set)
        D_vv_est = np.mean(D_vv_set)
        D_vv_err = np.std(D_vv_set, ddof=1) / np.sqrt(len(D_vv_set))
        print "D_vv = {} +/- {}".format(D_vv_est, D_vv_err)
    else:
        print "D_vv = {}".format(D_vv_set[0])

def ener_plot():
    E_t = np.loadtxt(args.enerfile)
    plt.plot(E_t[:,0], E_t[:,1])
    plt.plot(E_t[:,0], E_t[:,2])
    plt.legend(("Kinetic Energy", "Potential Energy"))
    plt.xlabel("Time, Tau")
    plt.title("Energy dynamics for rho={}".format(args.density))

if __name__ == "__main__":
    pos_plot()
    vel_plot()
    ener_plot()
