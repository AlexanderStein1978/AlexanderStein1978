#pragma once


class Oscillator
{
public:
    Oscillator();

    void initialize(const double omega, const double gamma, const double deltaT);

    double newValue(const double deltaAmplitude);

    double getAmplitudeFor(const double t) const;

    inline double getLastAmp() const
    {
        return mLastAmp;
    }

private:
    double mHOsq, mLastAmp, mLastV, mAmpAmp0Term, mAmpV0Term, mVAmp0Term, mVV0Term, mOmega, mGamma;
    static const double HALF;
};
