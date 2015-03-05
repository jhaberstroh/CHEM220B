import numpy as np
import argparse
import os
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
    
    if args.save is None:
        plt.show()
    else:
        outname = os.path.split(fname)[-1]
        outname = os.path.splitext(outname)[0]
        plt.savefig(args.save+"/"+outname+".png")
        plt.clf()
