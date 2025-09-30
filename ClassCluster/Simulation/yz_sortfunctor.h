//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef YZSORTFUNCTOR_H
#define YZSORTFUNCTOR_H


#include "particle.h"


class yz_SortFunctor
{
public:

    enum YorZ{yzY, yzZ};

    yz_SortFunctor(const Particle* const particles, const YorZ yz) : mParticles(particles), m_yz(yz)
    {
    }

    bool operator()(const int n1, const int n2)
    {
        if (n1 == -1) return false;
        if (n2 == -1) return true;
        return (m_yz == yzY ? mParticles[n1].R.Y() > mParticles[n2].R.Y() : mParticles[n1].R.Z() > mParticles[n2].R.Z());
    }

private:
    const Particle* const mParticles;
    const YorZ m_yz;
};

#endif // YZSORTFUNCTOR_H
