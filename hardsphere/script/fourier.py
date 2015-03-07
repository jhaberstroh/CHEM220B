import numpy as np
import argparse

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument("-file", default="fourier.csv", help="Path to csv files used in script")
parser.add_argument("-save", help="Path to save images out")
args = parser.parse_args()

if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

dat = np.loadtxt(args.file)
print dat.shape

n0 = 0
k = dat[n0:,0]
print k
coeff = dat[n0:,1:]

print coeff.shape
mean = np.array([np.mean(coeff[i,:]) for i in range(len(k))])
std  = np.array([np.std(coeff[i,:])  for i in range(len(k))])

N = 1000
plt.plot(k, mean*N)
plt.title("Mean of fourier distribution at different k")
plt.xlabel("k, units of 1/diameter")
plt.ylabel("Mean, Standard deviation")
plt.plot(k,  std*N)
plt.plot(k, -std*N)
safesaveplot(args.save,"fourier_mean_std.png")

plt.plot(k,  std**2*N)
plt.xlabel("k, units of 1/diameter")
plt.ylabel("Variance")
plt.title("Variance of fourier distribution at different k")
safesaveplot(args.save,"fourier_var.png")

values = [ np.pi, 2. * np. pi, 3. * np.pi ]
for val in values:
    k_ind = np.abs(k - val).argmin()
    density, bins = np.histogram(coeff[k_ind,:], bins=40, density=True)
    bin_ctr = (bins[:-1] + bins[1:]) / 2.
    print density.shape
    print bin_ctr.shape
    plt.plot(bin_ctr, np.log(density))
    plt.title("Fourier coefficient distribution for k = {}".format(val))
    plt.xlabel("k, units of 1/diameter")
    plt.ylabel("Log likelihood")

    mean = np.mean(coeff[k_ind,:])
    std  = np.std(coeff[k_ind,:])
    gauss = 1 / np.sqrt(2 * np.pi * std**2) * np.exp(-(bin_ctr - mean)**2 / (2 * std**2))
    plt.plot(bin_ctr, np.log(gauss))
    if args.save is None:
        plt.show()
    else:
        plt.savefig(args.save+"/fourier{}pi".format(int(round(val/np.pi)))+".png")
        plt.clf()
