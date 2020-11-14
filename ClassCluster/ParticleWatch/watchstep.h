#ifndef WATCHSTEP_H
#define WATCHSTEP_H


class ParticleContribution;


class WatchStep
{
public:
    WatchStep(const int iNumParticles);
    ~WatchStep();
    void reset();
    double get(const int particleIndex, const int coordinate) const;
    void set(const int particleIndex, const double x, const double y, const double z);
    void setSum(const double X, const double Y, const double Z);
    double getSumX() const;
    double getSumY() const;
    double getSumZ() const;

private:
    ParticleContribution* const contributions;
    const int numParticles;
    double sumX, sumY, sumZ;
};

#endif // WATCHSTEP_H
