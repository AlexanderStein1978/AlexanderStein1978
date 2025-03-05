#pragma once


class Oscillator
{
public:
    Oscillator();

    void initialize(const double omega, const double gamma, const double deltaT);

    double newValue(const double deltaAmplitude);

private:
    double mHOsq, mLastAmp, mLastV, mAmpAmp0Term, mAmpV0Term, mVAmp0Term, mVV0Term;
    static const double HALF;
};
