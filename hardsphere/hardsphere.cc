#include <stdio.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "MersenneTwister.h"
#include "parser.h"
#include "particles.h"


inline bool step_mc(double * particles, int N, double L, 
        MTRand & rng, double step_size)
{
    int atom = rng.randInt(N-1);
    assert(atom < N);
    assert(atom >= 0);
    double *dr = new double[3];
    dr[0] = rng.rand() * step_size - step_size / 2.;
    dr[1] = rng.rand() * step_size - step_size / 2.;
    dr[2] = rng.rand() * step_size - step_size / 2.;
    
    bool accept = true;
    for (int j = 0 ; j < N && accept ; j++)
    {
        double d = dist_ij(particles, atom, j, L, dr);
        if (d < 1.0 && atom != j)
        {
            accept=false;
        }
    }
    if (accept)
    {
        for (int dim = 0 ; dim < 3; dim ++)
        {
            particles[3*atom + dim] += dr[dim];
            particles[3*atom + dim] -= floor(particles[3*atom + dim] / L) * L;
        }
    }
    delete[] dr;
    return accept;
}

int main(int argc, char * argv[])
{
    int arg_i;
    int nsteq = 500; 
    int nstmc = 10000;
    int seed = 90210;
    double step_size = .3;
    double probe_rad = .45;
    int maxfouriernum = 1;
    double density = .5;
    int nstxyz = 10;
    int nstfourier = 10;
    // Parsing functions included from module parser.h, a small header 
    //    library that I wrote to make parsing command line arguments
    //    more transparent and modular.
    // Unfortunately, it does not automatically build a help file.
    for (arg_i = 1 ; arg_i < argc ; arg_i++)
    {
        parse_int(argc, argv, arg_i, "-nsteq", &nsteq, 
            "Number of sweeps for equilibration process");
        parse_int(argc, argv, arg_i, "-nstmc", &nstmc, 
            "Number of sweeps for production monte-carlo");
        parse_int(argc, argv, arg_i, "-seed", &seed, "Seed number");
        parse_double(argc, argv, arg_i, "-step_size", &step_size, 
            "Step size for MC moveset");
        parse_double(argc, argv, arg_i, "-density", &density, 
            "Density of simulation, in units of N/V");
        parse_double(argc, argv, arg_i, "-probe", &probe_rad, 
            "Radius for probe volume, units of particle diameters");
        #ifdef XYZ
        parse_int(argc, argv, arg_i, "-nstxyz", &nstxyz, 
            "Number of steps between performing position tracking");
        #endif //XYZ
        #ifdef FOURIER
        parse_int(argc, argv, arg_i, "-maxfouriernum", &maxfouriernum, 
            "Max fourier mode above the fundamental to record");
        parse_int(argc, argv, arg_i, "-nstfourier", &nstfourier, 
            "Max fourier mode above the fundamental to record");
        #endif //FOURIER
        help(argc, argv, arg_i);
    }

    int N_linear = 10;
    int N = N_linear * N_linear * N_linear;
    assert(density > 0.0 && density < 1.0);
    double L = double(N_linear) * cbrt(1./density);
      
    double* particles = new double [3 * N];
    initialize_grid_linarr(particles, N_linear, L);

    int nfouriervals = nstmc / nstfourier * 6;
    double* fouriers = new double[maxfouriernum * nfouriervals];
    
    #ifdef VERBOSE
    printf("nsteq: %d nstmc: %d\n", nsteq, nstmc);
    printf("L: %f\n", L);
    printf("math.pi: %f\n", M_PI);
    #endif
    
    #ifdef TEST
    if (!test_grid_linarr(particles, L, N))
    {
        printf("Linear grid test failed, halting execution\n");
        exit(1);
    }
    #endif

    int step;
    int acc;
    bool pass;
    MTRand rng(seed);

    // ============================================
    // Perform the equilibration Monte Carlo loop
    // ============================================
    acc = 0;
    // Steps are counted in units of number-of-sweeps out of convention
    for (step = 0 ; step < nsteq ; step ++)
    {   
        #ifdef ACCEPTANCE
        if (step % 50 == 0)
        {
            std::cout << "Step = " << step / N <<  " (" << step << ")" << std::endl;
            if (step > 0)
            {
                std::cout << "Acceptance rate: " <<  float(acc)/step/N * 100 << "%" << std::endl;
            }
        }
        #endif
        // step_mc is included step-by-step for modularity and clarity
        for (int step_i = 0 ; step_i < N ; step_i++)
        {
            pass = step_mc(particles, N, L, rng, step_size);
            if (pass)
            {
                acc++;
            }
        }

    }

    // ============================================
    // Perform the main Monte Carlo loop
    // ============================================
    acc = 0;
    for (step = 0 ; step < nstmc; step ++)
    {   
        for (int step_i = 0 ; step_i < N ; step_i++)
        {
            pass = step_mc(particles, N, L, rng, step_size);
            if (pass)
            {
                acc++;
            }
        }
        analyze_particles(particles, N, step, L,
            nstxyz, 0,                                  //xyz
            10, probe_rad,                              //probe volume
            nstfourier, fouriers, maxfouriernum, nfouriervals,  //fourier coeff
            50);                                        //g_r

        #ifdef ACCEPTANCE
        if (step % 50 == 0)
        {
            std::cout << "Step = " << step / N <<  " (" << step << ")" << std::endl;
            if (step > 0)
            {
                std::cout << "Acceptance rate: " <<  float(acc)/step / N * 100 << "%" << std::endl;
            }
        }
        #endif
    }

    // Post-processing
    #ifdef FOURIER
    postprocess_fourier(L, fouriers, maxfouriernum, nfouriervals);
    #endif
    delete[] fouriers;
    delete[] particles;
}
