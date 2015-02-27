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
        if bplot:
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
            if bplot:
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
            if bplot:
                print "Interact {}-{} = {}".format(i, j, interact_val)
            interact[i,j] = interact_val
            interact[j,i] = interact_val
    
    # Reintroduce units and include prefactors
    norm_R     = norm     * R_AW_A**3  
    overlap_R  = overlap  * R_AW_A**3 
    interact_R = interact / R_AW_A**3  * ( 4. * np.pi / (2 * np.pi)**3)
    return norm_R, overlap_R, interact_R

def compute_gr(sig_A_A = 2.7, bplot = False):
    sig_W_A = 2.7
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
    

    norm_R, overlap_R, interact_R = pc_coefficient_mtx(h_nl_interp, R_AW_A, 
                                                       r, basis_r, 
                                                       k, basis_k, 
                                                       bplot = bplot)

    if bplot:
        print norm_R
        print overlap_R
        print interact_R
    
    A = overlap_R + interact_R
    b = norm_R
    c_AW_i = LA.solve(A,-b)
    c_AW_r = np.zeros(N)
    c_AW_k = np.zeros(N)
    for i in xrange(4):
        c_AW_r += basis_r[i,:] * c_AW_i[i]
        c_AW_k += basis_k[i,:] * c_AW_i[i]


    h_AW_k = c_AW_k + (h_nl_interp * c_AW_k)

    dk_A = np.diff(k_A)
    N_gr = 400
    r_gr = np.linspace(-2, np.log10(R_AW_A + 10), N_gr)
    r_gr = np.power(10, r_gr)
    g_AW_r = np.interp(r_gr, r_A, c_AW_r)
    g_AW_r[r_gr > R_AW_A ] = 0
    g_AW_r2 = np.zeros(N_gr)
    h_AW_r  = np.zeros(N_gr)
    for i, r_i in enumerate(r_gr):
        integrand_k = h_nl_interp * c_AW_k * k_A * np.sin(k_A * r_i)
        g_AW_r2[i] = np.sum((integrand_k[1:] + integrand_k[:-1]) / 2. * dk_A) / (2 * np.pi **2 * r_i)

    
    h_AA_k = .033 * h_AW_k * c_AW_k
    h_AA_r = np.zeros(N_gr)
    for i, r_i in enumerate(r_gr):
        integrand_k = h_AA_k * k_A * np.sin(k_A * r_i)
        h_AA_r[i] = np.sum((integrand_k[1:] + integrand_k[:-1]) / 2. * dk_A) / (2 * np.pi **2 * r_i)
    h_AA_r[r_gr < sig_A_A] = 0

    return r_gr, g_AW_r + g_AW_r2 + 1, h_AA_r + 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot", action="store_true")
    args = parser.parse_args()
    
    l = []
    for sig_A in np.linspace(3, 5, 3):
        r, g_AW_r, g_AA_r = compute_gr(sig_A)
        R = (sig_A + 2.7) / 2.
        plt.plot(r - R, g_AW_r)
        l.append("sig={}".format(sig_A))
    plt.title("g_AW(r) at different radii")
    plt.legend(l)
    plt.show()
    
    N_R = 40
    R_max = 10.0
    sig_A_list = np.linspace(0, (R_max - 2.7), N_R)
    R_list = (sig_A_list + 2.7) / 2.0
    g_AW = np.zeros(N_R)
    for i, sig_i in enumerate(sig_A_list):
        r, gr, _ = compute_gr(sig_i)
        subset = (r > R_list[i])
        gr = gr[subset]
        r  = r [subset]
        r_min = r[0]
        # Use the magnitude of the slope
        g_slope = abs( (gr[1] - gr[0])/(r[1] - r[0]) )
        g_peak = gr[0] +  (r[0] - R_list[i]) * g_slope
        g_AW[i] = g_peak
    integrand = g_AW * 4 * np.pi * np.square(R_list)
    factor = .596 * .033
    delta_mu = np.cumsum((integrand[1:] + integrand[:-1]) / 2. * np.diff(R_list))
    delta_mu *= factor
    delta_mu = np.insert(delta_mu, 0, 0.)
    plt.plot(sig_A_list, delta_mu)
    plt.title("Solvation free energy for hydrophobe from Pratt-Chandler")
    plt.ylabel("Delta mu, KCal/mol")
    plt.xlabel("Solute diameter, Angstrom")
    plt.show()

    l = []
    for sig_A in np.linspace(3, 5, 3):
        r, g_AW_r, g_AA_r = compute_gr(sig_A)
        R = (sig_A + 2.7) / 2.
        plt.plot(r - sig_A, g_AA_r)
        l.append("sig={}".format(sig_A))
    plt.title("g_AA(r) at different radii")
    plt.legend(l)
    plt.show()

        
    
    

    
    


if __name__ == "__main__":
    main()
