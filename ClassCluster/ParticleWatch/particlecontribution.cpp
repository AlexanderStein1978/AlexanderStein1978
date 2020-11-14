#include "particlecontribution.h"

#include <cassert>


ParticleContribution::ParticleContribution() : x(0.0), y(0.0), z(0.0)
{
}

double ParticleContribution::get(const int coordinate) const
{
    switch (coordinate)
    {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            return 0.0;
            break;
    }
}

void ParticleContribution::reset()
{
    x = y = z = 0.0;
}

void ParticleContribution::set(const double ix, const double iy, const double iz)
{
    assert(x == 0.0 && y == 0.0 && z == 0.0);
    x = ix;
    y = iy;
    z = iz;
}
