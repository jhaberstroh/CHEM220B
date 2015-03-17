#include <iostream>
#include <cmath>
#include <cassert>
#include "MersenneTwister.h"
#include "parser.h"
#include "particles.h"
#include "velocities.h"


inline double forces_ij_lennardjones(double* particles, int pi, int pj, 
    double ljenergy, double ljdiameter, double ljcutoff, 
    double L, double* forces)
{
    double* dist_ij_dir = new double [3];
    double d = dist_ij(particles, pi, pj, L, NULL, dist_ij_dir);
    double f6 = pow(ljdiameter / d, 6);
    double fmag_ij = 0;
    double potential = 0;
    if (d < ljcutoff)
    {
        fmag_ij   = 4.*ljenergy * (-12.*f6*f6 + 6.*f6 ) / d;
        potential = 4.*ljenergy * (f6*f6 - f6);
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            forces[3 * pi + dim] -= fmag_ij * dist_ij_dir[dim];
            forces[3 * pj + dim] += fmag_ij * dist_ij_dir[dim];
        }
    }
    delete[] dist_ij_dir;
    return potential;
}

int step_mc_verlet(double* particles, double* velocity, int N, 
    double ljdiameter, double ljenergy, double ljcutoff, double dt, double L,
    double& potential_out)
{
    double potential = 0;
    double* particles_copy = new double [3*N];
    double* velocity_copy  = new double [3*N];
    double* forces = new double [3*N];
    for (int index = 0 ; index < 3 * N ; index++)
    {
        particles_copy[index] = particles[index];
        velocity_copy[index]  = velocity[index];
    }

    // Perform first half-step
    for (int index = 0 ; index < 3 * N ; index++)
    {
        forces[index]         = 0;
    }
    for (int pi = 0 ; pi < N ; pi++)
    {
        for (int pj = pi+1 ; pj < N ; pj++)
        {
            forces_ij_lennardjones(particles, pi, pj,  
                ljenergy, ljdiameter, ljcutoff, L, forces);
        }
        // Update half-velocities for particle pi
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            velocity_copy [3 * pi + dim] += forces[ 3 * pi + dim ] * dt / 2.;
            particles_copy[3 * pi + dim] += velocity_copy[ 3 * pi + dim ] * dt;
        }
    }
    
    // Check PBC
    for (int index = 0 ; index < 3 * N ; index++)
    {
        if (particles_copy[index] < -L)
        {
            std::cerr << "Particle far back at " << particles_copy[index] << std::endl;
        }
        if (particles_copy[index] > 2*L)
        {
            std::cerr << "Particle far up at " << particles_copy[index] 
              << " vs " << L << std::endl;
        }
        assert(particles_copy[index] >= -L && particles_copy[index] <= 2*L);
    }


    // Perform second half-step 
    //  -- changes marked with comments
    for (int index = 0 ; index < 3 * N ; index++)
    {
        forces[index]         = 0;
    }
    for (int pi = 0 ; pi < N ; pi++)
    {
        for (int pj = pi+1 ; pj < N ; pj++)
        {
            // Use particles_copy instead of particles
            potential += forces_ij_lennardjones(particles_copy, pi, pj, 
                ljenergy, ljdiameter, ljcutoff, L,  forces);
        }
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            // Only update velocities with t+dt forces
            velocity_copy [3 * pi + dim] += forces[ 3 * pi + dim ] * dt / 2.;
        }
    }

    // Copy results back into original arrays upon successful integration
    for (int index = 0 ; index < 3 * N ; index++)
    {
        particles[index] = fmod(particles_copy[index] + L, L);
        velocity [index] =  velocity_copy[index];
        assert(particles[index] >= 0 && particles[index] < L);
    }
    delete[] particles_copy;
    delete[] velocity_copy;
    delete[] forces;
    potential_out=potential;
}


int main(int argc, char * argv[])
{
    
    int nsteq = 500;
    int nstmd = 10000;
    int seed = 90210;
    double ljdiameter = 1;
    double ljenergy = 1.;
    double density = .8;
    double ljcutoff = 2.5;
    double dt = .001;
    int nstxyz = 10;
    double probe_rad = .45;
    int maxfouriernum = 1;
    double T = 3.0;

    for (int arg_i = 1 ; arg_i < argc ; arg_i++)
    {
        parse_int(argc, argv, arg_i, "-nsteq", &nsteq);
        parse_int(argc, argv, arg_i, "-nstmd", &nstmd);
        parse_int(argc, argv, arg_i, "-seed", &seed);
        parse_double(argc, argv, arg_i, "-ljdiameter", &ljdiameter);
        parse_double(argc, argv, arg_i, "-ljenergy", &ljenergy);
        parse_double(argc, argv, arg_i, "-ljcutoff", &ljcutoff);
        parse_double(argc, argv, arg_i, "-density", &density);
        parse_double(argc, argv, arg_i, "-dt", &dt);
        parse_int(argc, argv, arg_i, "-nstxyz", &nstxyz);
    }
    int nfouriervals = nstmd / 10 * 6;
    double* fouriers = new double[maxfouriernum * nfouriervals];
    
    MTRand rng(seed);
    int N_linear = 6;
    int N = N_linear * N_linear * N_linear;
    assert(density > 0.0 && density < 1.0);
    double L = double(N_linear) * cbrt(1./density);
      
    double* particles = new double [3 * N];
    double* velocity  = new double [3 * N];
    double Utot = 0;
    initialize_grid_linarr(particles, N_linear, L);
    generate_veldist(velocity, N, T, rng);
    remove_com(velocity, N);

    #ifdef VERBOSE
    std::cout << "Box size = " << L << ", ljcutoff = " << ljcutoff 
      << ", N = " << N << std::endl;
    #endif

    int step = 0;
    for (step = 0 ; step < nsteq ; step++)
    {
        #ifdef VERBOSE
        if (step % 50 == 0)
        {
            std::cout << "Step " << step << std::endl;
        }
        #endif

        double Utot_out = 0;
        double Ktot;
        step_mc_verlet(particles, velocity, N, ljdiameter, ljenergy, 
           ljcutoff, dt, L, Utot_out);
        Utot = Utot_out;
        Ktot = kinetic_energy(velocity, N);
        std::cerr << "Potential Energy: " << Utot;
        std::cerr << " Kinetic Energy: " << Ktot << std::endl;

        if (step % 50 == 0 && step != 0)
        {
            remove_com(velocity, N);
            thermostat_vrescale(velocity, N, T);
        }

        analyze_particles(particles, N, step, L,
            nstxyz, N,                               //xyz
            10, probe_rad,                              //probe volume
            10, fouriers, maxfouriernum, nfouriervals,  //fourier coeff
            50);                                        //g_r
    }

    delete[] fouriers;

    #ifdef VERBOSE
    std::cout << "Success!" << std::endl;
    #endif

}
