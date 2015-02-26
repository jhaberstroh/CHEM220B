import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-file", default="fourier.csv", help="Path to csv files used in script")
args = parser.parse_args()

dat = np.loadtxt(args.file)
print dat.shape

n0 = 0
k = dat[n0:,0]
print k
coeff = dat[n0:,1:]

print coeff.shape
mean = [np.mean(coeff[i,:]) for i in range(len(k))]
std  = [np.std(coeff[i,:])  for i in range(len(k))]

#plt.plot(2*np.pi/k, mean)
#plt.show()

plt.plot(k, std)
plt.show()

values = [ np.pi, 2. * np. pi, 3. * np.pi ]
for val in values:
    k_ind = np.abs(k - val).argmin()
    density, bins = np.histogram(coeff[k_ind,:], bins=40, density=True)
    bin_ctr = (bins[:-1] + bins[1:]) / 2.
    print density.shape
    print bin_ctr.shape
    plt.plot(bin_ctr, np.log(density))

    mean = np.mean(coeff[k_ind,:])
    std  = np.std(coeff[k_ind,:])
    gauss = 1 / np.sqrt(2 * np.pi * std**2) * np.exp(-(bin_ctr - mean)**2 / (2 * std**2))
    plt.plot(bin_ctr, np.log(gauss))
    plt.show()
