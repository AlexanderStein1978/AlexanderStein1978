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
    
    enum boundParticles{BoundAL = 6, NCandidates = 10};

    Particle() : next(nullptr), prev(nullptr), xp(0), yp(0), zp(0), NB(0), NC(0), MNB(0), WallPosIndex(-1), Fixed(false), WaveParticle(false) {}

    Particle *next, *prev;
    Binding bound[BoundAL], candidates[NCandidates];
    int xp, yp, zp, NB, NC, MNB, WallPosIndex;
    bool Fixed, WaveParticle;
    Vector R, v, lR, lv;
};

#endif // PARTICLE_H
