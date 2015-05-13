import numpy as np
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
parser.add_argument("-load", required=True, help="txt file to load")
parser.add_argument("-save", help="Location to save figures out")
args = parser.parse_args()

if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

data = np.loadtxt(args.load)
T_set = set(data[:,0])
dt_set = set(data[:,1])

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.0, 3.0)

T = .4
dt = .001
sub = (data[:, 0] == T) * (data[:, 1] == dt) * (data[:,2] < 200.)
plt.plot(data[sub, 2], data[sub, 3])
plt.xlabel("Time, reduced units")
plt.ylabel("Order parameter")
plt.title("Langevin timeseries, T=.4, dt=.1")
safesaveplot(args.save, "langevin_trajectory")

dt = .001    
for T in T_set:
    sub = (data[:, 0] == T) * (data[:, 1] == dt)
    min_prob = 1E-4
    bins = np.linspace(-1.7, 1.7, 100)
    x = bins[:-1] / 2. + bins[1:] / 2.
    y  = np.exp(-np.square(x+1)*np.square(x-1) / T)
    y /= np.sum( (x[1:]-x[:-1]) * (y[:-1]+ y[1:])*.5)
    bin_sub = (y > min_prob)
    x = x[bin_sub]
    y = y[bin_sub]
    hist, bins = np.histogram(data[sub, 3], bins = bins, density=True)
    plt.semilogy(x, hist[bin_sub])
    plt.semilogy(x, y)
    plt.ylim([1E-4, 1E1])
    plt.xlabel("Order parameter")
    plt.ylabel("Probability, log-scale")
    plt.title("Distribtuion, T={}, dt={}".format(T, dt))
    safesaveplot(args.save, "langevin_dist_{}T_{}dt".format(T, dt))

T = .2
for dt in dt_set:
    sub = (data[:, 0] == T) * (data[:, 1] == dt)
    min_prob = 1E-4
    bins = np.linspace(-1.7, 1.7, 100)
    x = bins[:-1] / 2. + bins[1:] / 2.
    y  = np.exp(-np.square(x+1)*np.square(x-1) / T)
    y /= np.sum( (x[1:]-x[:-1]) * (y[:-1]+ y[1:])*.5)
    bin_sub = (y > min_prob)
    x = x[bin_sub]
    y = y[bin_sub]
    hist, bins = np.histogram(data[sub, 3], bins = bins, density=True)
    plt.semilogy(x, hist[bin_sub])
    plt.semilogy(x, y)
    plt.ylim([1E-4, 1E1])
    plt.xlabel("Order parameter")
    plt.ylabel("Probability, log-scale")
    plt.title("Distribtuion, T={}, dt={}".format(T, dt))
    safesaveplot(args.save, "langevin_dist_{}T_{}dt".format(T, dt))
