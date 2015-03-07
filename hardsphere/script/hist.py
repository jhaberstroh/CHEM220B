import numpy as np
import argparse
import os

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()

parser = argparse.ArgumentParser()
parser.add_argument("-dir", default="", help="Path to csv files used in script")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()
if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

files=[args.dir + "/" + "small.csv",
       args.dir + "/" + "large-2.0.csv",
       args.dir + "/" + "large-3.0.csv",
       args.dir + "/" + "large-4.0.csv"]

fit=[False,
        True,
        True,
        True]
diameter=np.array([.9, 2, 3, 4])
means   =np.zeros(diameter.size)
stds    =np.zeros(diameter.size)

i=0
for fname, bFit in zip(files,fit):
    dat = np.loadtxt(fname)
    plt.title(fname)

    # HISTOGRAM
    # Create bins for all integers (and only integers) 
    MAX = int(max(dat))
    bins = np.linspace(-.5, MAX + .5, MAX + 2)
    hist,bins = np.histogram(dat, bins=bins, density=False)
    # Truncate statistically poor values
    hist[ hist < 100 ]  = 0
    # Cast integers to float to allow for normalization
    hist = hist.astype(np.float)
    means[i] = np.mean(dat)
    stds[i]  = np.std(dat)
    i += 1

    if not bFit:
        width = .1
        plt.bar(bins[:-1] - bins[0] - width/2., hist/np.sum(hist), width=width)
    
    if bFit:
        free_energy = np.log(hist / np.sum(hist))
        min_fe = np.log(100/np.sum(hist))
        plt.plot(free_energy, 'o')
        # FIT FUNCTION
        mean = np.mean(dat)
        std  = np.std(dat)
        N = 1. / np.sqrt(2. * np.pi * std**2)
        print "MEAN = {}".format(np.mean(dat))
        MIN = int(min(dat))
        xfit = np.linspace(MIN, MAX, 1000)
        yfit = np.exp(-np.square( xfit - mean ) / (2 * std **2) ) * N
        plt.plot(xfit, np.log(yfit))
        plt.vlines(mean, min_fe, 0)
        plt.ylim([min_fe, 0])
    
    outname = os.path.split(fname)[-1]
    outname = os.path.splitext(outname)[0]
    safesaveplot(args.save, outname+".png")

xmax = 4.5
x = np.linspace(0, xmax, 100)
r = x / 2.
density = .5
pv = 4./3. * np.pi * np.power(r,3) * density
plt.plot(x, pv)

plt.plot(diameter, means, 'o')
plt.vlines(diameter, means-stds, means+stds)
plt.xlim([0, xmax])
plt.title("Mean with diameter, compared to estimate")
plt.xlabel("Diameter, units of hardsphere radii")
safesaveplot(args.save, "/density-mean-var.png")

plt.plot(diameter, stds**2, 'o')
plt.title("Variance with diameter")
plt.xlabel("Diameter, units of hardsphere radii")
plt.xlim([0, xmax])
safesaveplot(args.save, "/density-var.png")
