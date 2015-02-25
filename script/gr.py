import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-file", default="fourier.csv", help="Path to csv files used in script")
args = parser.parse_args()

dat = np.loadtxt(args.file)
print dat.shape

density, bins = np.histogram(dat[dat != 0], bins=100)
bin_ctr = (bins[:-1] + bins[1:]) / 2.
plt.plot(bin_ctr, density / (4 * np.pi * np.square(bin_ctr)))
plt.xlim([0,max(bin_ctr)])
plt.savefig("/home/jhaberstroh/Dropbox/Physics/gr")
