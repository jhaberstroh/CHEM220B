import numpy as np
import argparse
import os

# Header for good configuration
def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="+", 
        help="Path to csv files used in script")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()
if args.save:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt


# Begin actual script:
#   Load data, plot RDF scaled down by 4pi R^2
print args.file
for fname in args.file:
    dat = np.loadtxt(fname)
    print dat.shape
    
    bins = np.linspace(1,5,41)
    density, bins = np.histogram(dat[dat != 0], bins=bins)
    bin_ctr = (bins[:-1] + bins[1:]) / 2.
    plt.plot(bin_ctr, density / (4 * np.pi * np.square(bin_ctr)))
    plt.xlim([0,max(bin_ctr)])
    
    basename=os.path.basename(fname)
    safesaveplot(args.save, basename + ".png")


