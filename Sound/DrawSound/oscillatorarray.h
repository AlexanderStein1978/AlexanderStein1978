//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once


class Oscillator;

class OscillatorArray
{
public:

    struct Results
    {
        double *time, *frequency, **data;
        int numTimeSteps;

        Results();
        Results(Results& right);
        Results& operator=(Results& right);
    };

    static const int NumOscillators;

    OscillatorArray(const int numTimeSteps, const double deltaT);
    ~OscillatorArray();

    void setNewValue(const int timeIndex, const double time, const double amplitude);

    inline Results& getResults()
    {
        return mResults;
    }

private:
    Oscillator* mOscillators;
    Results mResults;
    double mLastAmplitude;
};
