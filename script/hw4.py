import numpy as np
import numpy.linalg as LA
import numpy.fft as FFT
import matplotlib.pyplot as plt
import argparse




def pc_coefficient_mtx(h_nl_interp, R_AW_A, r, basis_r, k, basis_k, bplot = False ):
    # All basis functions asymptote, so set the smallest point to f(0) = f(kmin) after avoiding numerical failure
    dr = np.diff(r)
    dk = np.diff(k)
    if bplot:
        plt.plot(r, basis_r[0,:])
        plt.plot(r, basis_r[1,:])
        plt.plot(r, basis_r[2,:])
        plt.plot(r, basis_r[3,:])
        plt.show()
        plt.plot(k, basis_k[0,:])
        plt.plot(k, basis_k[1,:])
        plt.plot(k, basis_k[2,:])
        plt.plot(k, basis_k[3,:])
        plt.show()
    
    # ====================================================
    # Compute overlaps (no units)
    # ====================================================
    
    norm = np.zeros(4)
    for i in xrange(4):
        integrand_r = basis_r * 4. * np.pi * np.power(r,2)
        norm_val = np.sum((integrand_r[i,1:] + integrand_r[i,:-1]) / 2. * dr )
        print "Norm {} = {}".format(i, norm_val)
        norm[i] = norm_val
    
    # ====================================================
    # Compute overlaps (no units)
    # ====================================================
    
    overlap = np.zeros((4,4))
    for i in xrange(4):
        for j in xrange(i,4):
            overlap_fun_r = basis_r[i,:] * basis_r[j,:] * 4 * np.pi * np.power(r,2)
            if bplot:
                plt.plot(r,overlap_fun_r)
                plt.show()
            overlap_val = np.sum((overlap_fun_r[1:] + overlap_fun_r[:-1]) / 2. * dr)
            print "Overlap {}-{} = {}".format(i, j, overlap_val)
            overlap[i,j] = overlap_val
            overlap[j,i] = overlap_val
    
    # ====================================================
    # Compute interactions (accounting for units)
    # ====================================================
    
    
    interact = np.zeros((4,4))
    for i in xrange(4):
        for j in xrange(i,4):
            integrand_k = basis_k[i,:] * basis_k[j,:] * h_nl_interp * np.power(k,2)
            if bplot:
                plt.plot(k, integrand_k)
                plt.show()
            interact_val = np.sum((integrand_k[1:] + integrand_k[:-1]) / 2. * dk)
            print "Interact {}-{} = {}".format(i, j, interact_val)
            interact[i,j] = interact_val
            interact[j,i] = interact_val
    
    norm_R     = norm     * R_AW_A**3
    overlap_R  = overlap  * R_AW_A**3
    interact_R = interact / R_AW_A**3 * 4. / (2 * np.pi)**3

    return norm_R, overlap_R, interact_R

def compute_gr(sig_W_A = 2.7, bplot = False):
    sig_A_A = 2.7
    R_AW_A  = (sig_A_A + sig_W_A) / 2.

    # Set up basis
    N = 1001
    basis_k=np.zeros((4,N))
    basis_r=np.zeros((4,N))
    # k is in reduced units
    k = np.linspace(-2, 2.0, N)
    k = np.power(10,k)
    r = np.linspace(-7, 0, N)
    r = np.power(10,r)
    for n in xrange(4):
        basis_r[n,:] = np.power(r - 1, n)
    basis_k[0,:] = (np.sin(k) / np.power(k,3))       - ( np.cos(k) / np.power(k,2))
    basis_k[1,:] = (np.sin(k) / np.power(k,3))       + ( 2. * (np.cos(k) - 1) / np.power(k,4))
    basis_k[2,:] = (-6. * np.sin(k) / np.power(k,5)) + (( (2. * np.cos(k)) + 4.) / np.power(k,4))
    basis_k[3,:] = (-6. * np.sin(k) / np.power(k,5)) + ( 24. * (1 - np.cos(k)) / np.power(k,6)) - (6. / np.power(k,4)) 
    basis_k *= 4 * np.pi * R_AW_A**3
    k[0] = 0
    r[0] = 0
    
    r_A = r * R_AW_A
    k_A = k / R_AW_A

    # Load and interpolate narten-levy data
    k_nl = []
    h_nl = []
    with open("narten_levy.txt") as f:
        i = 0
        for l in f:
            if i >= 2:
                k_nl.append( float(l.split()[0]) )
                h_nl.append( float(l.split()[1]) )
            i += 1
        h_nl = np.array(h_nl)
        k_nl = np.array(k_nl)

    # Introduce a factor of R_AW to undimensionalize h_nl
    h_nl_interp = np.interp(k_A, k_nl, h_nl)
    if bplot:
        plt.plot(k_nl, h_nl)
        plt.plot(k, h_nl_interp)
        plt.xlim([min(k_nl), max(k_nl) + 5])
        plt.show()
    

    norm_R, overlap_R, interact_R = pc_coefficient_mtx(h_nl_interp, R_AW_A, r, basis_r, k, basis_k, bplot = bplot)

    print norm_R
    print overlap_R
    print interact_R
    
    A = overlap_R + interact_R
    b = norm_R
    c_AW_i = LA.solve(A,-b)
    print c_AW_i
    c_AW_r = np.zeros(N)
    c_AW_k = np.zeros(N)
    for i in xrange(4):
        c_AW_r += basis_r[i,:] * c_AW_i[i]
        c_AW_k += basis_k[i,:] * c_AW_i[i]
    #plt.plot(r_A, c_AW_r)
    #plt.title("C(r)")
    #plt.show()

    # Reintroduce units

    #plt.plot(r, c_AW_r)
    #plt.title("c_AW(r)")
    #plt.show()
    #plt.plot(k, c_AW_k)
    #plt.title("c_AW(k)")
    #plt.show()

    #g_r = c_AW_r
    #plt.plot(r, g_r)
    #plt.title("c(r) - [0,1] interpolatoin")
    #plt.show()
    #g_r2 = np.zeros(N)
    #for i, r_i in enumerate(r[1:]):
    #    integrand = h_nl_interp * c_AW_k * k * np.sin(k * r_i)
    #    g_r2[i+1] = np.sum((integrand[1:] + integrand[:-1]) / 2. * dk) / (2 * np.pi * r_i ** 2)
    #plt.plot(r, g_r)
    #plt.plot(r, g_r2)
    #plt.title("g(r) - [0,1] interpolation")
    #plt.ylim([-1, 10])
    #plt.show()

    #plt.plot(r_A, c_AW_r)
    #plt.title("c(r)")
    #plt.show()
    #plt.plot(k_A, c_AW_k)
    #plt.title("c(k)")
    #plt.show()

##  Code below has arbitrary factors included for "correctness", 
##  A safety pig is provided for your benefit 
## _._ _..._ .-',     _.._(`))
##'-. `     '  /-._.-'    ',/
##   )         \            '.
##  / _    _    |             \
## |  a    a    /              |
## \   .-.                     ;  
##  '-('' ).-'       ,'       ;
##     '-;           |      .'
##        \           \    /
##        | 7  .__  _.-\   \
##        | |  |  ``/  /`  /
##       /,_|  |   /,_/   /
##          /,_/      '`-'
##
    dk_A = np.diff(k_A)
    N_gr = 400
    r_gr = np.linspace(-2, np.log10(R_AW_A + 8), N_gr)
    r_gr = np.power(10, r_gr)
    g_r = np.interp(r_gr, r_A, c_AW_r)
    g_r[r_gr > R_AW_A ] = 0
    g_r2 = np.zeros(N_gr)
    for i, r_i in enumerate(r_gr):
        integrand_k = h_nl_interp * c_AW_k * k_A * np.sin(k_A * r_i) / np.pi # Include mystery pi???
        g_r2[i] = np.sum((integrand_k[1:] + integrand_k[:-1]) / 2. * dk_A) / (2 * np.pi **2 * r_i)

    return r_gr, g_r + g_r2 + 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot", action="store_true")
    args = parser.parse_args()
    
    r, gr = compute_gr(3.0)
    plt.plot(r, gr)
    plt.title("g(r)")
    plt.show()
    r, gr = compute_gr(4.0)
    plt.plot(r, gr)
    plt.title("g(r)")
    plt.show()
    r, gr = compute_gr(5.0)
    plt.plot(r, gr)
    plt.title("g(r)")
    plt.show()
    

    
    


if __name__ == "__main__":
    main()
