import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-dir", default="", 
        help="Path to csv files used in script")
args = parser.parse_args()

txt = np.loadtxt(args.dir + "/" + "atoms.xyz", delimiter=',')
print txt
for i in xrange(3):
    plt.plot(txt[:,i])
    plt.show()
