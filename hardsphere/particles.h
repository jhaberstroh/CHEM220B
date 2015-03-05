#ifndef __PARTICLES_H__
#define __PARTICLES_H__
#include <cmath>
#include <stdio.h>

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
    double period = 0, double* disp_i = NULL, double* dist_dir=NULL)
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
    if (dist_dir)
    {
        for (int dim = 0 ; dim < 3 ; dim++)
        {
            dist_dir[dim] = dr[dim]/dist;
        }
    }

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

void analyze_particles(double* particles, int N, int step, double L,
            int nstxyz, int particle_num,               //xyz
            int nstprobe, double probe_rad,             //probe volume
            int nstfourier, double* fouriers, 
            int maxfouriernum, int nfouriervals,        //fourier coeff
            int nstgr)                                  //g_r
{

    #ifdef XYZOUT
    if (step % nstxyz == 0)
    {
        printf("%f, %f, %f\n", particles[3 * particle_num], 
                               particles[3 * particle_num + 1], 
                               particles[3 * particle_num + 2]);
    }
    #endif
    #ifdef SMALLSPHERE
    if (step % nstprobe == 0)
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
    if (step % nstprobe == 0)
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
    int offset = (step / nstfourier) * 6;
    if ((step % nstfourier == 0) && 
        (offset < nfouriervals))
    {
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
    if (step % nstgr == 0)
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
}

void postprocess_fourier(double L, double* fouriers, int maxfouriernum, 
    int nfouriervals)
{
    for (int k = 0 ; k < maxfouriernum ; k++)
    {
        std::cout << fourier_coeff(k,L) << " ";
        for (int nk = 0 ; nk < nfouriervals ; nk++)
        {
            std::cout << fouriers[k * nfouriervals + nk] << " ";
        }
        std::cout << std::endl;
    }
}


#endif //__PARTICLES_H__
