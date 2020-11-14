#include "watchstep.h"
#include "particlecontribution.h"


WatchStep::WatchStep(const int iNumParticles) : contributions(new ParticleContribution[iNumParticles]), numParticles(iNumParticles), sumX(0.0), sumY(0.0), sumZ(0.0)
{
}

WatchStep::~WatchStep()
{
    delete[] contributions;
}

void WatchStep::reset()
{
    for (int i=0; i < numParticles; ++i) contributions[i].reset();
    sumX = sumY = sumZ = 0.0;
}

double WatchStep::get(const int particleIndex, const int coordinate) const
{
    return contributions[particleIndex].get(coordinate);
}

void WatchStep::set(const int particleIndex, const double x, const double y, const double z)
{
    contributions[particleIndex].set(x, y, z);
}

void WatchStep::setSum(const double X, const double Y, const double Z)
{
    sumX = X;
    sumY = Y;
    sumZ = Z;
}

double WatchStep::getSumX() const
{
    return sumX;
}

double WatchStep::getSumY() const
{
    return sumY;
}

double WatchStep::getSumZ() const
{
    return sumZ;
}
