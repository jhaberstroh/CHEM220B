import numpy as np
import argparse

def safesaveplot(savedir=None, name=None):
    if savedir is None:
        plt.show()
    else:
        plt.savefig(savedir+'/'+name)
        plt.clf()

parser = argparse.ArgumentParser()
parser.add_argument("-save", help="Number of K points to use")
args = parser.parse_args()

if not args.save is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt

w = np.linspace(-3, 3, 300)
gammaz = (.1, .5, 1.)
for gamma in gammaz:
    S = gamma / ((gamma**2 / 4.) + np.square(w + 1)) - \
        gamma / ((gamma**2 / 4.) + np.square(w - 1)) 
    plt.plot(w, S)
plt.legend(["gamma = {}".format(x) for x in gammaz])
plt.xlabel("omega")
plt.ylabel("Response Function")
plt.title("Response function for various gammas, Langevin equation")
safesaveplot(args.save, "FP_response.png")

