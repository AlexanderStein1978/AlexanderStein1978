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
        double l1(mParticles[n1].v.lengthSquared()), l2(mParticles[n2].v.lengthSquared());
        if (l1 < l2 && l1 > l2) return false;
        return (l1 > l2);
    }

private:
    const Particle* const mParticles;
};

#endif // VSORTFUNCTOR_H
