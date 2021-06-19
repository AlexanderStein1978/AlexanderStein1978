#include "watchpoint.h"
#include "watchstep.h"


WatchPoint::WatchPoint(const int iNumSteps, const int iNumParticles) : steps(new WatchStep*[iNumSteps]), numSteps(iNumSteps)
{
    for (int i=0; i < numSteps; ++i) steps[i] = new WatchStep(iNumParticles);
}

WatchPoint::~WatchPoint()
{
    for (int i=0; i < numSteps; ++i) delete steps[i];
    delete[] steps;
}

void WatchPoint::reset()
{
    for (int i=0; i < numSteps; ++i) steps[i]->reset();
}

double WatchPoint::get(const int step, const int particleIndex, const int coordinate) const
{
    return steps[step]->get(particleIndex, coordinate);
}

void WatchPoint::set(const int step, const int particleIndex, const Vector &r)
{
    steps[step]->set(particleIndex, r);
}

void WatchPoint::setSum(const int step, const Vector &R)
{
    steps[step]->setSum(R);
}

double WatchPoint::getSumX(const int step) const
{
    return steps[step]->getSumX();
}

double WatchPoint::getSumY(const int step) const
{
    return steps[step]->getSumY();
}

double WatchPoint::getSumZ(const int step) const
{
    return steps[step]->getSumZ();
}
