#include "oscillator.h"

#include <cmath>


const double Oscillator::HALF = 0.5;


Oscillator::Oscillator() : mHOsq(0.0), mLastAmp(0.0), mLastV(0.0), mAmpAmp0Term(0.0), mAmpV0Term(0.0), mVAmp0Term(0.0), mVV0Term(0.0), mOmega(0.0), mGamma(0.0)
{
}

void Oscillator::initialize(const double omega, const double gamma, const double deltaT)
{
    const double expGammaT = exp(-1.0 * gamma * deltaT), coswt = cos(omega * deltaT), sinwt = sin(omega * deltaT), dOmega = 1.0 / omega, gammaDOmega = gamma * dOmega;
    mHOsq = HALF * omega * omega;
    mAmpAmp0Term = expGammaT * (coswt + gammaDOmega * sinwt);
    mAmpV0Term = expGammaT * dOmega * sinwt;
    mVAmp0Term = -1.0 * expGammaT * (omega + gamma * gammaDOmega) * sinwt;
    mVV0Term = expGammaT * (coswt - gammaDOmega * sinwt);
    mOmega = omega;
    mGamma = gamma;
}

double Oscillator::newValue(const double deltaAmplitude)
{
    const double currAmp = mLastAmp * mAmpAmp0Term + mLastV * mAmpV0Term;
    mLastV = mLastAmp * mVAmp0Term + mLastV * mVV0Term;
    mLastAmp = currAmp - deltaAmplitude;
    return HALF * mLastV * mLastV + mHOsq * mLastAmp * mLastAmp;
}

double Oscillator::getAmplitudeFor(const double t) const
{
    return exp(-1.0 * mGamma * t) * (mLastAmp * cos(mOmega * t) + (mLastAmp * mGamma + mLastV) * sin(mOmega * t) / mOmega);
}
