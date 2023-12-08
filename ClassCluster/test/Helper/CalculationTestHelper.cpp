#include "CalculationTestHelper.h"
#include "Calculation.h"
#include "particle.h"

#include <cmath>


CalculationTestHelper::CalculationTestHelper(Calculation* calc) : mCalc(calc)
{
}

bool CalculationTestHelper::bindsParticleTo(Particle* P1, Particle* P2) const
{
    for (int n=0; n < P1->NB; ++n) if (P1->bound[n].p == P2) return true;
    return false;
}

bool CalculationTestHelper::areParticlesBound(Particle* P1, Particle* P2) const
{
    return (bindsParticleTo(P1, P2) && bindsParticleTo(P2, P1));
}

bool CalculationTestHelper::areParticlesBound(const int index1, const int index2) const
{
    return areParticlesBound(mCalc->P + index1, mCalc->P + index2);
}

bool CalculationTestHelper::areParticleBindingsCorrectlyInitialized(int index) const
{
    return (areParticlesBound(mCalc->P + index, mCalc->P + index - 2)
         && areParticlesBound(mCalc->P + index, mCalc->P + index + 2)
         && areParticlesBound(mCalc->P + index, mCalc->P + index - 40)
         && areParticlesBound(mCalc->P + index, mCalc->P + index - 38)
         && areParticlesBound(mCalc->P + index, mCalc->P + index + 38)
         && areParticlesBound(mCalc->P + index, mCalc->P + index + 40));
}

void CalculationTestHelper::swapParticlePositions(const int index1, const int index2)
{
    std::swap(mCalc->P[index1].R, mCalc->P[index2].R);
}

void CalculationTestHelper::updateBindings()
{
    mCalc->UpdateBindings();
}

int CalculationTestHelper::getPXS() const
{
    return mCalc->PXS;
}

int CalculationTestHelper::getPZS() const
{
    return mCalc->PZS;
}

double CalculationTestHelper::getDist(const int index1, const int index2) const
{
    return Calculation::dist(mCalc->P + index1, mCalc->P + index2);
}

int CalculationTestHelper::getBoundParticleIndex(const int particleIndex, const int bindingIndex) const
{
    return mCalc->P[particleIndex].bound[bindingIndex].p - mCalc->P;
}

int CalculationTestHelper::getNumParticles() const
{
    return mCalc->N;
}

bool CalculationTestHelper::isNoBindingDoubled(const int particleIndex) const
{
    return !(mCalc->isBindingDoubled(particleIndex));
}

int CalculationTestHelper::getParticleNB(const int particleIndex) const
{
    return mCalc->P[particleIndex].NB;
}

int CalculationTestHelper::getParticleMNB(const int particleIndex) const
{
    return mCalc->P[particleIndex].MNB;
}

double CalculationTestHelper::getRandom()
{
    return static_cast<double>(rand()) / RAND_MAX;
}

double CalculationTestHelper::createRandomSpeedDistribution()
{
    double energy = 0.0;
    for (int n=0; n < mCalc->N; ++n)
    {
        mCalc->P[n].v.setX(0.5 - getRandom());
        mCalc->P[n].v.setY(0.5 - getRandom());
        mCalc->P[n].v.setZ(0.5 - getRandom());
        mCalc->P[n].v *= (1e4 * getRandom());
        energy += mCalc->P[n].v.lengthSquared();
    }
    return 0.5 * energy;
}

void CalculationTestHelper::setParticles(const int N, const Vector* const R, const Vector* const v, const int* const MNB)
{
    for (int x=0; x < mCalc->XS; x++) for (int y=0; y < mCalc->YS; y++) for (int z=0; z < mCalc->ZS; ++z) mCalc->G[x][y][z] = nullptr;
    mCalc->N = N;
    mCalc->mRandPOF1 = static_cast<double>(N) / (static_cast<double>(RAND_MAX) + 1.0);
    mCalc->mRandPOF2 = static_cast<double>(N-1) / (static_cast<double>(RAND_MAX) + 1.0);
    Vector F(double(mCalc->XS) / mCalc->MaxX, double(mCalc->YS) / mCalc->MaxY, double(mCalc->ZS) / mCalc->MaxZ * 1.001);
    for (int n=0; n<N; ++n)
    {
        mCalc->initializeParticle(mCalc->P[n], 5, 5, R[n], F, true);
        mCalc->P[n].v = v[n];
        mCalc->P[n].MNB = MNB[n];
    }
}

void CalculationTestHelper::addParticleBinding(const int n1, const int n2)
{
    mCalc->P[n1].bound[mCalc->P[n1].NB].p = mCalc->P + n2;
    mCalc->P[n2].bound[mCalc->P[n2].NB].p = mCalc->P + n1;
    mCalc->P[n1].bound[mCalc->P[n1].NB++].lastDist = mCalc->P[n2].bound[mCalc->P[n2].NB++].lastDist = (mCalc->P[n1].R - mCalc->P[n2].R).length();
}

Vector CalculationTestHelper::getCenterOfMass()
{
    Vector C;
    for (int n=0; n < mCalc->N; ++n) C += mCalc->P[n].R;
    return (C / mCalc->N);
}

Vector CalculationTestHelper::getAverageV()
{
    Vector v;
    for (int n=0; n < mCalc->N; ++n) v += mCalc->P[n].v;
    return (v / mCalc->N);
}

Vector CalculationTestHelper::getAngularMomentum(const Vector& C)
{
    Vector L;
    for (int n=0; n < mCalc->N; ++n) L += (mCalc->P[n].R - C).cross(mCalc->P[n].v);
    return L;
}

double CalculationTestHelper::getBindingAngle(const int leftIndex, const int centerIndex, const int rightIndex) const
{
    const Vector leftDiff = mCalc->P[leftIndex].R - mCalc->P[centerIndex].R, rightDiff = mCalc->P[rightIndex].R - mCalc->P[centerIndex].R;
    return acos(leftDiff.dot(rightDiff) / (leftDiff.length() * rightDiff.length()));
}

double CalculationTestHelper::getSpeedSum() const
{
    double result = 0.0;
    for (int n=0; n < mCalc->N; ++n) result += mCalc->P[n].v.length();
    return result;
}

void CalculationTestHelper::run(const int maxIteration)
{
    mCalc->mMaxIt = maxIteration;
    mCalc->start();
    mCalc->wait();
}

bool CalculationTestHelper::getU(const int n1, const int n2, double& U, const Vector *const t0, int pos, Vector* a, const bool collectCandidates) const
{
    return Calculation::Success == mCalc->getU(mCalc->P + n1, mCalc->P + n2, U, t0, static_cast<Calculation::Positions>(pos), a, collectCandidates);
}
