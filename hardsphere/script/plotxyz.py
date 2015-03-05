import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-dir", default="", 
        help="Path to csv files used in script")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()
if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

txt = np.loadtxt(args.dir + "/" + "atoms.xyz", delimiter=',')
print txt
for i in xrange(3):
    T = txt.shape[0]
    plt.plot(txt[:T/5,i])
    if args.save is None:
        plt.show()
    else:
        plt.savefig(args.save+"/"+"xyz{}.png".format(i))
        plt.clf()
