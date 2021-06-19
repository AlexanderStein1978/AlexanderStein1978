#ifndef PARTICLECONTRIBUTION_H
#define PARTICLECONTRIBUTION_H


#include "vector.h"


class ParticleContribution
{
public:
    ParticleContribution();

    void reset();
    void set(const Vector &i);
    double get(const int coordinate) const;

private:
    Vector r;
};

#endif // PARTICLECONTRIBUTION_H
