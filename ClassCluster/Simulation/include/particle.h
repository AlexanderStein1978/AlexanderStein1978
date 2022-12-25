#ifndef PARTICLE_H
#define PARTICLE_H


#include "vector.h"

struct Particle
{
    struct Binding
    {
        Binding() : p(nullptr), lastDist(0.0) {}

        Particle* p;
        double lastDist;
    };
    
    enum boundParticles{BoundAL = 4, NCandidates = 10};

    Particle() : next(nullptr), prev(nullptr), xp(0), yp(0), zp(0), NB(0), NC(0), MNB(0), E(0.0), deltaE(0.0), T(0.0), deltaT(0.0), U(0.0), deltaU(0.0) {}

    Particle *next, *prev;
    Binding bound[BoundAL], candidates[NCandidates];
    int xp, yp, zp, NB, NC, MNB;
    double E, deltaE, T, deltaT, U, deltaU;
    Vector R, v, lR, lv, aa, corr;
};

#endif // PARTICLE_H
