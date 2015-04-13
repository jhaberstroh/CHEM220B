from __future__ import division
from __future__ import print_function
import numpy as np
import numpy.random as rand
import scipy.spatial.distance as distance
import argparse

def xyzout(particles, filename, writemode='a', comment=None):
    N = particles.shape[0]
    with open(filename, writemode) as f:
        print(N, file=f)
        if comment is None:
            print("Comment blank", file=f)
        else:
            print(comment, file=f)
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
    parser = argparse.ArgumentParser()
    parser.add_argument("-nstmd", default=500, type=int, 
            metavar="STEPS",
            help="Number of steps to run for production MD.")
    parser.add_argument("-nsteq", default=500, type=int,
            metavar="STEPS",
            help="Numer of steps to run for equilibration & field MD.")
    parser.add_argument("-field", nargs=2, type=float, default=(0,0),
            metavar=("N", "STRENGTH"),
            help='''Pair of args for applied field, (1) coefficient of pi/L 
            for field harmonic (2) field strength.''') 
    parser.add_argument("-nstxyz", default=10, type=int,
            metavar="STEPS",
            help="Frequency with which to write to posfile.")
    parser.add_argument("-seed", default=90210, type=int,
            metavar="INT",
            help="RNG seed value.")
    parser.add_argument("-posfile", required=True,
            metavar="NAME",
            help="File to write XYZ data to.")
    parser.add_argument("-enerfile", required=True,
            metavar="NAME",
            help="File to write energy data to.")
    parser.add_argument("-density", required=True, type=float,
            metavar="FRACTION",
            help="Number density per cubic diameter")
    parser.add_argument("-N_linear", required=True, type=int,
            metavar="N",
            help="Number of particles per linear dimension")
    parser.set_defaults(nstvrescale=100, 
                        nstener=100,
                        nstvel=0,
                        kT=1.5)
    args = parser.parse_args()
    if args.density >= 1.0 or args.density < 0:
        raise ValueError("-density must be a value between 0.0 and 1.0")

    nsteq = args.nsteq
    nstmd = args.nstmd
    nstvrescale = args.nstvrescale
    nstxyz = args.nstxyz
    nstvel = args.nstvel
    nstener = args.nstener
    enerwrite = 'w'

    U_t = np.zeros(nstener)
    K_t = np.zeros(nstener)

    density = args.density
    N_linear = args.N_linear
    N        = N_linear ** 3
    rand.seed(args.seed)
    kT       = args.kT
    dt       = .001
    posfile  = args.posfile
    enerfile = args.enerfile
    
    L = 1. / np.power(density, 1/3) * N_linear
    L_minforce = np.power(2, 1/6) * N_linear
    opt_meshsize = min(L_minforce, L)
    particles = np.zeros((N,3))
    one_d_mesh = np.linspace(0, opt_meshsize, N_linear, endpoint=False)
    x, y, z = np.meshgrid( one_d_mesh, one_d_mesh, one_d_mesh )
    particles[:,0] = x.flatten()
    particles[:,1] = y.flatten()
    particles[:,2] = z.flatten()
    k_field = args.field[0] * np.pi / L
    f_field = args.field[1]

    velocities = rand.normal( 0, kT/2., (N,3))
    velocities -= np.mean(velocities, axis=0)
    for step in xrange(nsteq):
        if step == 0:
            F_ir, _ = ljforces(particles, L)
            F_ir[:, 0] += f_field * np.sin(k_field * particles[:, 0]) * k_field / np.sqrt(N)
            xyzout(particles, posfile, writemode='w', comment="t={} L={} k={} F={}".format(
                    dt * (step - nsteq), L, k_field, f_field))

        velocities += F_ir * dt/2
        particles  += velocities * dt
        pbc(particles, L)
        F_ir, U = ljforces(particles, L)
        F_ir[:, 0] += f_field * np.sin(k_field * particles[:, 0]) * k_field / np.sqrt(N)
        U += np.sum(f_field * np.cos(k_field * particles[:, 0])) / np.sqrt(N)
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
                scale = np.sqrt(3 * kT * (N-1)  / msv)
            print("Rescaling: {}".format(scale))
            velocities *= scale
        if nstxyz > 0 and step % nstxyz == 0 and step != 0:
            xyzout(particles, posfile, comment="t={} L={} k={} F={}".format(
                    dt * (step - nsteq), L, k_field, f_field))

    for step in xrange(nstmd):
        if step == 0:
            F_ir, _ = ljforces(particles, L)
            xyzout(particles, posfile, comment="t={} L={} k={} F={}".format(
                    dt * step, L, k_field, 0))

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
            xyzout(particles, posfile, comment="t={} L={} k={} F={}".format(
                    dt * step, L, k_field, 0))
        if nstvel > 0 and step % nstvel == 0 and step != 0:
            raise NotImplementedError("No velocity output has been implemented")


if __name__ == "__main__":
    main()
