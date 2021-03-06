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
parser.add_argument("-N", default=1000, type=int, help="Number of particles")
parser.add_argument("-T", default=100, type=int, help="Number of frames")
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
    outname = os.path.split(fname)[-1]
    outname = os.path.splitext(outname)[0]
    ind = outname.index("_")
    density = float("." + outname[ind+1:])
    
    N_frames = args.T
    # Half particles used becasue factor of 2 for symmetry in distance has
    #  already been accounted for
    half_particles = args.N/2

    bins = np.linspace(1,5,101)
    N_R, bins = np.histogram(dat[dat != 0], bins=bins)
    N_R = N_R.astype(float)
    N_R /= float(N_frames * half_particles)
    bin_ctr  = (bins[1:] + bins[:-1]) / 2.
    bin_size = (bins[1:] - bins[:-1])
    plt.plot(bin_ctr, N_R / (4 * np.pi * np.square(bin_ctr) * bin_size * density) )
    plt.xlim([0,max(bin_ctr)])
    plt.ylim([0,4.5])
    plt.ylabel("g(r)")
    plt.xlabel("r, particle diameters")
    plt.title("Density = {}".format(density))
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(3.0, 2.5)
    
    outname = os.path.split(fname)[-1]
    outname = os.path.splitext(outname)[0]
    safesaveplot(args.save, outname + ".png")


