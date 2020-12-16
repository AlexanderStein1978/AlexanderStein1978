#ifndef VSORTFUNCTOR_H
#define VSORTFUNCTOR_H


#include "particle.h"


class VSortFunctor
{
public:
    explicit VSortFunctor(const Particle* const particles) : mParticles(particles)
    {
    }

    bool operator()(const int n1, const int n2)
    {
        if (n1 == -1) return false;
        if (n2 == -1) return true;
        return  mParticles[n1].vX * mParticles[n1].vX + mParticles[n1].vY * mParticles[n1].vY + mParticles[n1].vZ * mParticles[n1].vZ >
                mParticles[n2].vX * mParticles[n2].vX + mParticles[n2].vY * mParticles[n2].vY + mParticles[n2].vZ * mParticles[n2].vZ;
    }

private:
    const Particle* const mParticles;
};

#endif // VSORTFUNCTOR_H
