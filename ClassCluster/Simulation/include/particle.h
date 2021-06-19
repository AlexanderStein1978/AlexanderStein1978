#ifndef PARTICLE_H
#define PARTICLE_H


#include "vector.h"

struct Particle
{
    enum boundParticles{NBound = 4};
    Particle *next, *prev, *bound[NBound];
    int xp, yp, zp, NB;
    double E, deltaE, T, deltaT, U, deltaU;
    Vector R, v, lR, lv, aa, corr;
};

#endif // PARTICLE_H
