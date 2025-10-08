//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef WATCHSTEP_H
#define WATCHSTEP_H


#include "vector.h"


class ParticleContribution;


class WatchStep
{
public:
    WatchStep(const int iNumParticles);
    ~WatchStep();
    void reset();
    double get(const int particleIndex, const int coordinate) const;
    void set(const int particleIndex, const Vector& r);
    void setSum(const Vector& sum);
    double getSumX() const;
    double getSumY() const;
    double getSumZ() const;

private:
    ParticleContribution* const contributions;
    const int numParticles;
    Vector sum;
};

#endif // WATCHSTEP_H
