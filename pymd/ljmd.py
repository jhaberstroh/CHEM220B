from __future__ import division
from __future__ import print_function
import numpy as np
import numpy.random as rand
import scipy.spatial.distance as distance


def xyzout(particles, filename, writemode='a'):
    N = particles.shape[0]
    with open(filename, writemode) as f:
        print(N, file=f)
        print("Comment blank", file=f)
        vdwsig = 3.4
        for i in xrange(particles.shape[0]):
            print("C     " + 
                    "{:>8.5f} ".format(particles[i,0] * vdwsig) +
                    "{:>8.5f} ".format(particles[i,1] * vdwsig) + 
                    "{:>8.5f}". format(particles[i,2] * vdwsig), file=f)

def enerout(U_t, K_t, step, dt, nstener, filename, writemode='a'):
    with open(filename, writemode) as f:
        for t, (U_i,K_i) in enumerate(zip(U_t,K_t)):
            T =  dt * (step - (nstener-1) + t)
            print("{:>10.4f} ".format(T),
                  "{:>10.4f} ".format(U_i),
                  "{:>10.4f} ".format(K_i), file=f)


def pbc(particles, L):
    particles[particles >= L] -= L
    particles[particles < 0 ] += L

def ljforces(particles, L):
    N = particles.shape[0]
    # Compute all distance vectors from particle i to particle j
    #   with a minimum image convention
    Dx_ijr   = particles[np.newaxis, :, :] - particles[:, np.newaxis, :]
    Dx_ijr[Dx_ijr  >  L/2] -= L
    Dx_ijr[Dx_ijr <= -L/2] += L

    # Compute several powers of the distance magnitudes
    D2_ij  = np.sum(np.square(Dx_ijr), axis = 2)
    D_ij  =  np.reciprocal(np.sqrt(D2_ij))
    D6_ij  = np.power(D_ij, 6)
    D12_ij = np.square(D6_ij)

    # Normalize all distance vectors and set the diagonals to zero
    Dx_ijr /= np.sqrt(D2_ij[:, :, np.newaxis])
    np.fill_diagonal(Dx_ijr[:,:,0], 0)
    np.fill_diagonal(Dx_ijr[:,:,1], 0)
    np.fill_diagonal(Dx_ijr[:,:,2], 0)
    
    # Compute forces for all pairs of particles within range of one
    #   another, and zero all others. Also zero self-interactions
    Fmag_ij = np.zeros( (N, N) )
    U_ij    = np.zeros( (N, N) )
    inrange = (D2_ij < 2.5**2)  *  (D2_ij != 0)
    Fmag_ij[inrange] = ( - 48 * D12_ij[inrange]   \
                         + 24 *  D6_ij[inrange] ) \
                         * D_ij[inrange]
    U_ij[inrange] =   4 *  ( +  D12_ij[inrange]   \
                             -   D6_ij[inrange] ) 
    # Point forces along the distance vector direction and sum
    F_ijr  = Dx_ijr * Fmag_ij[:, :, np.newaxis]
    F_ir = np.sum(F_ijr, axis = 1)
    return F_ir, .5 * np.sum(U_ij, axis=(0,1)) 

def main():
    nsteq = 300
    nstmd = 500
    nstvrescale = 100
    nstxyz = 10
    nstvel = 0
    nstener = 100
    enerwrite = 'w'

    U_t = np.zeros(nstener)
    K_t = np.zeros(nstener)

    density = .3
    N_linear = 6
    N        = N_linear ** 3
    rand.seed(90210)
    kT       = 1.5
    dt       = .001
    posfile  = "mini_lj.xyz"
    enerfile = "mini_lj.ener"
    
    L = 1. / np.power(density, 1/3) * N_linear
    L_minforce = np.power(2, 1/6) * N_linear
    opt_meshsize = min(L_minforce, L)
    particles = np.zeros((N,3))
    one_d_mesh = np.linspace(0, opt_meshsize, N_linear, endpoint=False)
    x, y, z = np.meshgrid( one_d_mesh, one_d_mesh, one_d_mesh )
    particles[:,0] = x.flatten()
    particles[:,1] = y.flatten()
    particles[:,2] = z.flatten()

    velocities = rand.normal( 0, kT/2., (N,3))
    velocities -= np.mean(velocities, axis=0)
    for step in xrange(nsteq):
        if step == 0:
            F_ir, _ = ljforces(particles, L)
            xyzout(particles, posfile, writemode='w')

        velocities += F_ir * dt/2
        particles  += velocities * dt
        pbc(particles, L)
        F_ir, U = ljforces(particles, L)
        velocities += F_ir * dt/2

        U_t[step % nstener] = U
        K_t[step % nstener] = .5 * np.sum(np.square(velocities), axis=(0,1))
        if step % nstener == nstener - 1:
            enerout(U_t, K_t, step-nsteq, dt, nstener, enerfile, writemode=enerwrite)
            enerwrite='a'
        if step % nstvrescale == 0:
            msv = np.sum(np.square(velocities), axis=(0,1))
            if msv == 0:
                scale = 0
            else:
                scale = np.sqrt(3/2 * kT * (N-1)  / msv)
            print("Rescaling: {}".format(scale))
            velocities *= scale
        if nstxyz > 0 and step % nstxyz == 0 and step != 0:
            xyzout(particles, posfile)

    for step in xrange(nstmd):
        if step == 0:
            F_ir, _ = ljforces(particles, L)
            xyzout(particles, posfile)

        velocities += F_ir * dt/2
        particles += velocities * dt
        pbc(particles, L)
        F_ir, U = ljforces(particles, L)
        U_t[step % nstener] = U
        K_t[step % nstener] = .5 * np.sum(np.square(velocities), axis=(0,1))
        velocities += F_ir * dt/2

        if step % nstener == nstener - 1:
            enerout(U_t, K_t, step, dt, nstener, enerfile, writemode=enerwrite)
            enerwrite='a'
        if nstxyz > 0 and step % nstxyz == 0 and step != 0:
            xyzout(particles, posfile)
        if nstvel > 0 and step % nstvel == 0 and step != 0:
            xyzout(particles, posfile)

if __name__ == "__main__":
    main()
