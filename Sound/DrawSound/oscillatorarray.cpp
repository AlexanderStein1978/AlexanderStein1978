#include "oscillatorarray.h"
#include "oscillator.h"
#include "utils.h"


namespace
{
    const double gamma = 0.1;
    const double MinOmega = 10;
    const double MaxOmega = 5000;
    const double OmegaStep = (MaxOmega - MinOmega) / (NumOscillators - 1);
}

const int OscillatorArray::NumOscillators = 100;


OscillatorArray::OscillatorArray(const int numTimeSteps) : mOscillators(new Oscillator[NumOscillators])
{
    mResults.time = new double[numTimeSteps];
    mResults.frequency = new double[NumOscillators];
    mResults.data = Create(numTimeSteps, NumOscillators);
    mResults.numTimeSteps = numTimeSteps;

    int n;
    double omega;
    for (n=0, omega = MinOmega; n < NumOscillators; ++n, omega += OmegaStep)
    {
        mResults.frequency[n] = omega;
        mOscillators[n].initialize(omega, gamma);
    }
}

OscillatorArray::~OscillatorArray()
{
    delete[] mOscillators;
    if (nullptr != mResults.data)
    {
        delete[] mResults.time;
        delete[] mResults.frequency;
        Destroy(mResults.data, numTimeSteps);
    }
}

void OscillatorArray::setNewValue(const int timeIndex, const double time, const double amplitude)
{
    mResults.time[timeIndex] = time;
    for (int n=0; n < NumOscillators; ++n) mResults.data[timeIndex][n] = mOscillators[n].newValue(time, amplitude);
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
    data = nullptr;
    return *this;
}
