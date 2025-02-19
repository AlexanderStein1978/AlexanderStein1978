#include "oscillator.h"

#include <cmath>


double Oscillator::HALF = 0.5;
double Oscillator::M1 = -1.0;


Oscillator::Oscillator() : mOmega(0.0), mGamma(0.0), mLastAmplitude(0.0), mLastOA(0.0), mLastV(0.0), mGmO(0.0),
    mGpO(0.0), md2O(0.0), mLastAvA(0.0), mhOsq(0.0), mLastExpGmO(1.0), mLastExpGpO(1.0)
{
}

void Oscillator::initialize(const double omega, const double gamma)
{
    mOmega = omega;
    mGamma = gamma;
    mGmO = gamma - omega;
    mGpO = gamma + omega;
    md2O = HALF / omega;
    mhOsq = HALF * omega * omega;
}

double Oscillator::newValue(const double amplitude, const double time)
{
    double avA(HALF * (amplitude + mLastAmplitude)), deltaAvA(avA - mLastAvA), oA = mLastOA - deltaAvA;
    double c1 = (oA + (oA * mGmO + mLastV) * md2O) / mLastExpGmO;
    double c2 = M1 * (oA * (mGmO + mLastV) * md2O) / mLastExpGpO;
    mLastAmplitude = amplitude;
    mLastAvA = avA;
    mLastExpGmO = exp(M1 * mGmO * time);
    mLastExpGpO = exp(M1 * mGpO * time);
    mLastOA = c1 * mLastExpGmO + c2 * mLastExpGpO;
    mLastV = M1 * (c1 * mGmO * mLastExpGmO + c2 * mGpO * mLastExpGpO);
    return HALF * mLastV * mLastV + mhOsq * mLastOA * mLastOA;
}
