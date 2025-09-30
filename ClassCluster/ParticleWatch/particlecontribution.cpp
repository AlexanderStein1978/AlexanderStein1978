//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "particlecontribution.h"

#include <cassert>


ParticleContribution::ParticleContribution() : r()
{
}

double ParticleContribution::get(const int coordinate) const
{
    switch (coordinate)
    {
        case 0:
            return r.X();
        case 1:
            return r.Y();
        case 2:
            return r.Z();
        default:
            return 0.0;
            break;
    }
}

void ParticleContribution::reset()
{
    r = Vector();
}

void ParticleContribution::set(const Vector &i)
{
    assert(r.X() == 0.0 && r.Y() == 0.0 && r.Z() == 0.0);
    r = i;
}
