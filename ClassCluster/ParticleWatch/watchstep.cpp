//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "watchstep.h"
#include "particlecontribution.h"


WatchStep::WatchStep(const int iNumParticles) : contributions(new ParticleContribution[iNumParticles]), numParticles(iNumParticles), sum()
{
}

WatchStep::~WatchStep()
{
    delete[] contributions;
}

void WatchStep::reset()
{
    for (int i=0; i < numParticles; ++i) contributions[i].reset();
    sum = Vector();
}

double WatchStep::get(const int particleIndex, const int coordinate) const
{
    return contributions[particleIndex].get(coordinate);
}

void WatchStep::set(const int particleIndex, const Vector &r)
{
    contributions[particleIndex].set(r);
}

void WatchStep::setSum(const Vector& Sum)
{
    sum = Sum;
}

double WatchStep::getSumX() const
{
    return sum.X();
}

double WatchStep::getSumY() const
{
    return sum.Y();
}

double WatchStep::getSumZ() const
{
    return sum.Z();
}
