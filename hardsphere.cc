#include <stdio.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "MersenneTwister.h"
#include "parser.h"

double inline fourier_coeff(int k, double L)
{
    return 2.0 * double(k+1) * M_PI / L;
}

void initialize_grid_linarr(double* particles, int N_linear, double L)
{
    double space = L / double(N_linear);
    for (int i = 0 ; i < N_linear ; i++)
    {
        for (int j = 0 ; j < N_linear ; j++)
        {
            for (int k = 0 ; k < N_linear ; k++)
            {
                int ind = N_linear * N_linear * i +
                                     N_linear * j +
                                                k;
                particles[3 * ind    ] = space * (double) i;
                particles[3 * ind + 1] = space * (double) j;
                particles[3 * ind + 2] = space * (double) k;
            }
        }
    }
}

// Computes nearest replica between i and j with optional displacement 
//   of i 
inline double dist_ij(double* particles, int i, int j, 
    double period = 0, double* disp_i = NULL)
{
    double* dr = new double[3];
    for (int dim = 0 ; dim < 3 ; dim++)
    {
        dr[dim] = particles[(3*i) + dim] - particles[(3*j) + dim];
    }
    if (disp_i)
    {
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            dr[dim] += disp_i[dim];
        }
    }

    // Cube nearest image convention
    if (period)
    {
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            dr[dim] -= round( dr[dim] / period ) * period;
        }
    }

    double dist = 0;
    for (int dim = 0 ; dim < 3 ; dim++)
    {
        dist += dr[dim] * dr[dim];
    }
    dist = sqrt(dist);

    delete[] dr;
    return dist;
}

bool test_grid_linarr(double* particles, double L, int N)
{
    // Check that all particles replicated on top
    //   of one another with box size == 1 
    for (int i = 0 ; i < N ; i++)
    {
        if (dist_ij(particles, 0, i, L / cbrt(N) ) > .01)
        {
            printf("Fail 0 %f\n", L / cbrt(double(N)));
            return false;
        }
    }

    for (int i = 0 ; i < N ; i ++)
    {
        for (int j = i+1 ; j < N ; j++)
        {
            double dist = dist_ij(particles, i, j, L);
            if (dist < .9)
            {
                printf("Fail 1 %d %d %f\n", i, j, dist);
                return false;
            }
            if (dist > (2 * L))
            {
                printf("Fail 2\n");
                return false;
            }
        }
    }
    return true;
}

int probe_volume(double* particles, int N, double* probe_r, double probe_rad)
{
    double * dist = new double[3];
    
    int occ = 0;
    for (int j = 0 ; j < N ; j++)
    {
        double d_tot = 0;
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            dist[dim] = probe_r[dim] - particles[3*j + dim];
            d_tot += dist[dim] * dist[dim];
        }
        d_tot = sqrt(d_tot);
        if (d_tot < probe_rad)
        {
            occ++;
        }
    }
    delete[] dist;
    return occ;
}

inline bool step_mc(double * particles, int N, double L, 
        MTRand & rng, double step_size)
{
    int atom = rng.randInt(N-1);
    assert(atom < N);
    assert(atom >= 0);
    double *dr = new double[3];
    dr[0] = rng.rand() * step_size;
    dr[1] = rng.rand() * step_size;
    dr[2] = rng.rand() * step_size;
    
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
    for (arg_i = 1 ; arg_i < argc ; arg_i++)
    {
        parse_int(argc, argv, arg_i, "-nsteq", &nsteq);
        parse_int(argc, argv, arg_i, "-nstmc", &nstmc);
        parse_int(argc, argv, arg_i, "-seed", &seed);
        parse_double(argc, argv, arg_i, "-step_size", &step_size);
        parse_double(argc, argv, arg_i, "-probe", &probe_rad);
        parse_int(argc, argv, arg_i, "-maxfouriernum", &maxfouriernum);
        parse_double(argc, argv, arg_i, "-density", &density);
    }

    int N_linear = 10;
    int N = N_linear * N_linear * N_linear;
    assert(density > 0.0 && density < 1.0);
    double L = double(N_linear) * cbrt(1./density);
      
    double* particles = new double [3 * N];
    initialize_grid_linarr(particles, N_linear, L);

#ifdef FOURIER
    int nfouriervals = nstmc / 10 * 6;
    double* fouriers = new double[maxfouriernum * nfouriervals];
#endif
    
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

    // Perform the equilibration Monte Carlo loop
    acc = 0;
    for (step = 0 ; step < (nsteq * N); step ++)
    {   
#ifdef VERBOSE
        if (step % (50 * N) == 0)
        {
            std::cout << "Step = " << step / N <<  " (" << step << ")" << std::endl;
            if (step > 0)
            {
                std::cout << "Acceptance rate: " <<  float(acc)/step * 100 << "%" << std::endl;
            }
        }
#endif
        pass = step_mc(particles, N, L, rng, step_size);
        if (pass)
        {
            acc++;
        }

    }

    // Perform the main Monte Carlo loop
    acc = 0;
    for (step = 0 ; step < (nstmc * N); step ++)
    {   
        pass = step_mc(particles, N, L, rng, step_size);

#ifdef VERBOSE
        if (step % (50 * N) == 0)
        {
            std::cout << "Step = " << step / N <<  " (" << step << ")" << std::endl;
            if (step > 0)
            {
                std::cout << "Acceptance rate: " <<  float(acc)/step * 100 << "%" << std::endl;
            }
        }
#endif
#ifdef ACCEPTANCE
        if (step % (50 * N) == 0)
        {
            std::cout << "Step = " << step / N <<  " (" << step << ")" << std::endl;
            if (step > 0)
            {
                std::cout << "Acceptance rate: " <<  float(acc)/step * 100 << "%" << std::endl;
            }
        }
#endif
#ifdef XYZOUT
        if (step % (50 * N) == 0)
        {
            printf("%f, %f, %f\n", particles[3 * 555], 
                                    particles[3 * 555 + 1], 
                                    particles[3 * 555 + 2]);
        }
#endif
#ifdef SMALLSPHERE
        if (step % (10 * N) == 0)
        {
            double * probe_r = new double[3];
            probe_r[0] = 5.0;
            probe_r[1] = 5.0;
            probe_r[2] = 5.0;
            int occ = probe_volume(particles, N, probe_r, .45);
            printf("%d ",occ);
            delete[] probe_r;
        }

#endif

#ifdef LARGESPHERE
        if (step % (10 * N) == 0)
        {
            double * probe_r = new double[3];
            for (int i = 2 ; i < L - 2 ; i++)
            {
                for (int j = 2 ; j < L - 2 ; j++)
                {
                    for (int k = 2 ; k < L - 2 ; k++)
                    {

                        probe_r[0] = i;
                        probe_r[1] = j;
                        probe_r[2] = k;
                        int occ = probe_volume(particles, N, probe_r, probe_rad);
                        printf("%d ",occ);
                    }
                }
            }
            delete[] probe_r;
        }
#endif

#ifdef FOURIER
        if ((step % (10 * N) == 0) && 
            (step / 10 / N * 6 < nfouriervals))
        {
            int offset = (step / 10 / N) * 6;
            for (int k = 0 ; k < maxfouriernum ; k++)
            {
                double kval = fourier_coeff(k,L);
                int k_i = (k * nfouriervals) + offset;

                for (int i = 0 ; i < 6 ; i++)
                {
                    fouriers[k_i + i] = 0;
                }
                for (int i = 0 ; i < N ; i++)
                {
                    for (int j = 0 ; j < 3 ; j++)
                    {
                        fouriers[k_i+j  ] += cos(particles[3*i+j] * kval);
                        fouriers[k_i+j+3] += sin(particles[3*i+j] * kval);
                    }
                }
                for (int i = 0 ; i < 6 ; i++)
                {
                    fouriers[k_i + i] /= double(N);
                }
            }
        }
#endif

#ifdef GR
        if (step % (50 * N) == 0)
        {
            for (int i = 0 ; i < N ; i++)
            {
                for (int j = i+1 ; j < N ; j++)
                {
                    double gr_dist = dist_ij(particles, i, j, L);
                    if (gr_dist < L/2.0)
                    {
                        std::cout << gr_dist << " ";
                    }
                }
            }
        }

#endif

        if (pass)
        {
            acc++;
        }
    }

    // Post-processing

#ifdef FOURIER
    for (int k = 0 ; k < maxfouriernum ; k++)
    {
        std::cout << fourier_coeff(k,L) << ",";
        for (int nk = 0 ; nk < nfouriervals ; nk++)
        {
            std::cout << fouriers[k * nfouriervals + nk] << " ";
        }
        std::cout << std::endl;
    }
    delete[] fouriers;
#endif
    delete[] particles;
}
