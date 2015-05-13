import numpy as np
import scipy.special as spec
import argparse
import string

def safesaveplot(savedir=None, name=None, suffix='.png'):
    if not name is None:
        name = string.replace(name, '.', '_') + suffix
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()

parser = argparse.ArgumentParser()
parser.add_argument("-save", help="Location to save figures out")
parser.add_argument("-name", help="Special name for output figure")
parser.add_argument("-harmonic", action='store_true', help="Check the harmonic appx. Overrides -traj argument")
parser.add_argument("-traj", nargs="+", help="Trajectory-based splitting rate files")
parser.add_argument("-traj_q0", type=float, nargs="+", help="q0 values for trajectory-based splitting rate files")
parser.add_argument("-kT", type=float, required=True)
args = parser.parse_args()

if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

fig = plt.gcf()
fig.set_size_inches(6.0, 4.0)


kT = args.kT
q = np.linspace(-1, 1, 2001)
w = np.exp( (q+1)**2 * (q-1)**2 / kT)
integrand = (w[1:] + w[:-1])/2. * (q[1:] - q[:-1])
phi = np.cumsum(integrand) / np.sum(integrand)
plt.plot(q[1:], phi)
plt.title("Comparison of methods, kT={}".format(args.kT))
plt.xlabel("q")
plt.ylabel("Splitting probability")
plt.tight_layout()
if args.harmonic:
    harm = 1 - .5 * spec.erfc(np.sqrt(2 / kT) * q)
    plt.plot(q, harm)
    plt.legend(['Onsager result', 'harmonic result'], loc=2)
    safesaveplot(args.save, "harmonic_kT{}".format(kT))

if not args.traj is None:
    print args.traj
    phi = np.zeros(len(args.traj))
    for i, file in enumerate(args.traj):
        with open(file) as f:
            line = f.readline()
            B_counts = line.count("B")
            norm     = len(line)
            phi[i] = float(B_counts) / len(line)
    plt.plot(args.traj_q0, phi)
    plt.legend(['Onsager result', 'trajectory result'], loc=2)
    if args.name:
        safesaveplot(args.save, args.name)
    else:
        safesaveplot(args.save, "trajectories_kT{}".format(kT))
