import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="+", 
        help="Path to csv files used in script")
args = parser.parse_args()

print args.file
for fname in args.file:
    dat = np.loadtxt(fname)
    print dat.shape
    basename=os.path.basename(fname)
    
    bins = np.linspace(1,5,41)
    density, bins = np.histogram(dat[dat != 0], bins=bins)
    bin_ctr = (bins[:-1] + bins[1:]) / 2.
    plt.plot(bin_ctr, density / (4 * np.pi * np.square(bin_ctr)))
    plt.xlim([0,max(bin_ctr)])
    plt.savefig("/home/jhaberstroh/Dropbox/Physics/gr/" + basename + ".png")
    plt.clf()


