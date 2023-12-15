#pragma once


#include "vector.h"
#include "Calculation.h"


struct Particle;


class CalculationTestHelper
{
public:
    CalculationTestHelper(Calculation* calc);
    ~CalculationTestHelper();

    bool bindsParticleTo(Particle* P1, Particle* P2) const;
    bool areParticlesBound(Particle* P1, Particle* P2) const;
    bool areParticlesBound(const int index1, const int index2) const;
    bool areParticleBindingsCorrectlyInitialized(int index) const;
    void swapParticlePositions(const int index1, const int index2);
    void updateBindings();
    int getPXS() const;
    int getPZS() const;
    double getDist(const int index1, const int index2) const;
    int getBoundParticleIndex(const int particleIndex, const int bindingIndex) const;
    int getNumParticles() const;
    bool isNoBindingDoubled(const int particleIndex) const;
    int getParticleNB(const int particleIndex) const;
    int getParticleMNB(const int particleIndex) const;
    static double getRandom();
    double createRandomSpeedDistribution();
    void setParticles(const int N, const Vector* const R, const Vector* const v, const int* const MNB);
    void addParticleBinding(const int n1, const int n2);
    Vector getCenterOfMass();
    Vector getAverageV();
    Vector getAngularMomentum(const Vector& C);
    Vector getPositionDifference(const int index1, const int index2) const;
    void run(const int maxIteration);
    double getBindingAngle(const int leftIndex, const int centerIndex, const int rightIndex) const;
    double getSpeedSum() const;
    bool getU(const int n1, const int n2, double &U, const Vector* const t0, int pos, Vector *a, const bool collectCandidates) const;
    void replacePotential(const Calculation::PotRole role, double* const newPot, double* const newdPdR);
    void resetPotential(const Calculation::PotRole role);

    inline int getNumPotentialPoints()
    {
        return mCalc->NPot;
    }

private:
    Calculation* mCalc;
    double **Pot, **dPdR;
};
