import numpy as np
import numpy.fft as FFT
import argparse

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument("-velfile", required=True, help="Path to csv files used in script")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()
if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt


def main():
    vel = np.loadtxt(args.velfile)
    vel = vel.reshape((vel.shape[0], vel.shape[1]/3, 3))
    print vel
    plt.plot(vel[:,0,0])
    plt.show()
    Ct = np.zeros((vel.shape[0]/2))
    for i in xrange(vel.shape[1]):
        for j in xrange(vel.shape[2]):
            f_w = FFT.rfft(vel[:,i,j])
            mag = f_w * np.conj(f_w)
            Ct_ij = FFT.irfft(mag)
            Ct += Ct_ij[:len(Ct_ij)/2] / vel.shape[0] / np.pi * 2.
    plt.plot(Ct / vel.shape[1] / vel.shape[2])
    plt.show()




if __name__ == "__main__":
    main()
