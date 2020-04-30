#ifndef DELTAESORTFUNCTOR_H
#define DELTAESORTFUNCTOR_H


#include "particle.h"


class DeltaESortFunctor
{
public:
    explicit DeltaESortFunctor(const Particle* const particles) : mParticles(particles)
    {
    }

    bool operator()(const int n1, const int n2)
    {
        if (n1 == -1) return false;
        if (n2 == -1) return true;
        return mParticles[n1].deltaE < mParticles[n2].deltaE;
    }

private:
    const Particle* const mParticles;
};

#endif // DELTAESORTFUNCTOR_H
