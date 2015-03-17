#include <cmath>
#include <iostream>

int thermostat_vrescale(double* velocity, int N, double T)
{
    double vtot_sq = 0;
    for (int index = 0 ; index < 3 * N ; index++)
    {
        vtot_sq += velocity[index] * velocity[index];
    }
    // Perform the rescaling to the thermal ensemble
    double rescale;
    if (vtot_sq == 0)
    {
        rescale = 0;
    }
    else
    {
        rescale = sqrt(3. / 2. * T * (N-1) ) / sqrt(vtot_sq);
        std::cerr << "Rescale factor: " << rescale << std::endl;
    }
    for (int index = 0 ; index < 3 * N ; index++)
    {
        velocity[index] *= rescale;
    }
}

double kinetic_energy(double* velocity, int N)
{
    double Ktot = 0;
    for (int index = 0 ; index < 3 * N ; index++)
    {
        Ktot += velocity[index] * velocity[index];
    }
    return Ktot;
}


int generate_veldist(double* velocity, int N, double T, MTRand& rng)
{
    for (int index = 0 ; index < 3 * N ; index++)
    {
        velocity[index] = rng.rand(T/2.) - T/4.;
    }
}


int remove_com(double* velocity, int N)
{
    double* vtot = new double[3];
    for (int i = 0 ; i < 3 ; i ++)
    {
        vtot[i] = 0;
    }
    for (int index = 0 ; index < 3 * N ; index++)
    {
        vtot[index%3] += velocity[index];
    }
    for (int index = 0 ; index < 3 * N ; index++)
    {
        velocity[index] -= vtot[index%3] / N;
    }
    delete[] vtot;
}
