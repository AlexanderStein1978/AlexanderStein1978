#include "oscillatorarray.h"
#include "oscillator.h"
#include "utils.h"

#include <math.h>

const int OscillatorArray::NumOscillators = 200;

namespace
{
    const double Gamma = 2.0;
    const double MinF = 10;
    const double MaxF = 10000;
    const double FStep = (MaxF - MinF) / (OscillatorArray::NumOscillators - 1);
}

OscillatorArray::OscillatorArray(const int numTimeSteps, const double deltaT) : mOscillators(new Oscillator[NumOscillators])
{
    mResults.time = new double[numTimeSteps];
    mResults.frequency = new double[NumOscillators];
    mResults.data = Create(numTimeSteps, NumOscillators);
    mResults.numTimeSteps = numTimeSteps;

    int n;
    double f, timeStep = deltaT / (numTimeSteps - 1);
    for (n=0, f = MinF; n < NumOscillators; ++n, f += FStep)
    {
        mResults.frequency[n] = f;
        mOscillators[n].initialize(2.0 * M_PI * f, Gamma, timeStep);
    }
}

OscillatorArray::~OscillatorArray()
{
    delete[] mOscillators;
    if (nullptr != mResults.data)
    {
        delete[] mResults.time;
        delete[] mResults.frequency;
        Destroy(mResults.data, mResults.numTimeSteps);
    }
}

void OscillatorArray::setNewValue(const int timeIndex, const double time, const double amplitude)
{
    mResults.time[timeIndex] = time;
    double deltaAmp = amplitude - mLastAmplitude;
    mLastAmplitude = amplitude;
    for (int n=0; n < NumOscillators; ++n) mResults.data[timeIndex][n] = mOscillators[n].newValue(deltaAmp);
}


OscillatorArray::Results::Results() : time(nullptr), frequency(nullptr), data(nullptr), numTimeSteps(0)
{
}


OscillatorArray::Results::Results(Results& right)
{
    *this = right;
}

OscillatorArray::Results& OscillatorArray::Results::operator=(Results& right)
{
    time = right.time;
    right.time = nullptr;
    frequency = right.frequency;
    right.frequency = nullptr;
    data = right.data;
    right.data = nullptr;
    numTimeSteps = right.numTimeSteps;
    right.numTimeSteps = 0;
    return *this;
}
