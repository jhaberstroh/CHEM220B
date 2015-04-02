import numpy as np
import numpy.fft as FFT
import argparse
from StringIO import StringIO

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument("-velfile", required=True, help="Path to velocity data")
parser.add_argument("-xyzfile", required=True, help="Path to xyz data")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()
if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt


def loadxyz(filename):
    xyz = []
    with open(filename) as f:
        i = 0
        data = ""
        for l in f:
            if i==0:
                n_particle = int(l)
            elif i==1:
                comment = l
            else:
                data = data + l
            i+=1
            if i==(n_particle+2):
                i=0
                data = StringIO(data)
                frame = np.loadtxt(data, usecols=[1,2,3])
                xyz.append(frame)
                data = ""
    xyz = np.array(xyz)
    return xyz
            


def main():
    L = 6.4633
    xyz = loadxyz(args.xyzfile)
    print xyz.shape
    shift = L/2. - xyz[0,:,:] 
    xyz = (xyz + shift + L) % L
    xyz -= L/2.
    for i in xrange(xyz.shape[1]):
        plt.plot(xyz[:,i,0])
        plt.plot(xyz[:,i,1])
        plt.plot(xyz[:,i,2])
    safesaveplot(args.save, "x_ensemble.png")
    
    plt.plot(xyz[:,0,0])
    safesaveplot(args.save, "x_particle.png")
    exit()


    vel = np.loadtxt(args.velfile)
    vel = vel.reshape((vel.shape[0], vel.shape[1]/3, 3))
    print vel
    plt.plot(vel[:,0,0])
    safesaveplot(args.save, "v_particle.png")
    Ct = np.zeros((vel.shape[0]/2))
    for i in xrange(vel.shape[1]):
        for j in xrange(vel.shape[2]):
            f_w = FFT.rfft(vel[:,i,j])
            mag = f_w * np.conj(f_w)
            Ct_ij = FFT.irfft(mag)
            Ct += Ct_ij[:len(Ct_ij)/2] / vel.shape[0] * 2
    plt.plot(Ct / vel.shape[1] / vel.shape[2])
    safesaveplot(args.save, "v_autocorr.png")




if __name__ == "__main__":
    main()
