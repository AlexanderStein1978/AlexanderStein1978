#pragma once

class Oscillator
{
public:
    Oscillator();

    void initialize(const double omega, const double gamma);

    double newValue(const double amplitude, const double time);

private:
    double mOmega, mGamma, mLastAmplitude, mLastOA, mLastV, mGmO, mGpO, md2O, mLastAvA, mhOsq, mLastExpGmO, mLastExpGpO;
    static double HALF, M1;
};
