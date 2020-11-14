#ifndef PARTICLECONTRIBUTION_H
#define PARTICLECONTRIBUTION_H


class ParticleContribution
{
public:
    ParticleContribution();

    void reset();
    void set(const double ix, const double iy, const double iz);
    double get(const int coordinate) const;

private:
    double x, y, z;
};

#endif // PARTICLECONTRIBUTION_H
