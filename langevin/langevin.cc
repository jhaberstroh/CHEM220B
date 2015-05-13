#include "MersenneTwister.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "parser.h"

// From Wikipedia Box-Muller Page (2015-05-08 12:42:00)
double generateGaussian(MTRand& rng, double mu, double sigma)
{
  const double epsilon = std::numeric_limits<double>::min();
  const double two_pi = 2.0*3.14159265358979323846;
 
  static double z0, z1;
  static bool generate;
  generate = !generate;
 
  if (!generate)
     return z1 * sigma + mu;
 
  double u1, u2;
  do
   {
     u1 = rng.rand();
     u2 = rng.rand();
   }
  while ( u1 <= epsilon );
 
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}


int main(int argc, char** argv)
{
    double q0 = 0.0;
    double kT_split = .1;
    int seed = 90210;
    for (int arg_i = 1 ; arg_i < argc ; arg_i++)
    {
        parse_double(argc, argv, arg_i, "-q0", &q0, 
            "Starting position");
        parse_double(argc, argv, arg_i, "-kT_split", &kT_split, 
            "Temperature for splitting calculation");
        parse_int(argc, argv, arg_i, "-seed", &seed, 
            "Seed value for random numbers");
    }
    MTRand rng(seed);

    int nsteps = 1E7;
    double q_t = q0;

    #ifdef SPLITTING_RATE
    double dt = .001;
    double kT = kT_split;
    double D = 1.;
    double g  = kT / D;
    for (int step = 0 ; step < nsteps ; step++)
    {
        double F = 4. * q_t * (q_t + 1.) * (q_t - 1.);
        double R = generateGaussian(rng, 0, sqrt(2 * D * dt));
        double q_t1 = q_t - (1./g) * F * dt + R;
        q_t = q_t1;
        if (q_t >= 1)
        {
            std::cout << "B";
            exit(0);
        }
        if (q_t <= -1)
        {
            std::cout << "A";
            exit(0);
        }
    }
    
    #else
    const double dt_list[] = {.001, .002, .003};
    const double kT_list[] = {.4, .2, .1, .05};
    for (int dt_i = 0 ; dt_i < 3 ; dt_i++)
    {
        double dt = dt_list[dt_i];
        for (int kT_i = 0 ; kT_i < 4 ; kT_i++)
        {
            double kT = kT_list[kT_i];
            double D = 1.;
            double g  = kT / D;
            for (int step = 0 ; step < nsteps ; step++)
            {
                double F = 4. * q_t * (q_t + 1.) * (q_t - 1.);
                double R = generateGaussian(rng, 0, sqrt(2 * D * dt));
                double q_t1 = q_t - (1./g) * F * dt + R;
                q_t = q_t1;
                if (step % 100 == 0)
                {
                    std::cout << kT << " " << dt << " " << 
                            step*dt << " " << q_t << std::endl;
                }
            }
        }
    }
    #endif 
}
