#ifndef PARTICLE_H
#define PARTICLE_H


#include "vector.h"

struct Particle
{
    struct Binding
    {
        Particle* p;
        double lastDist;
    };
    
    enum boundParticles{BoundAL = 4, NCandidates = 10};
    Particle *next, *prev;
    Binding bound[BoundAL], candidates[NCandidates];
    int xp, yp, zp, NB, NC, MNB;
    double E, deltaE, T, deltaT, U, deltaU;
    Vector R, v, lR, lv, aa, corr;
};

#endif // PARTICLE_H
