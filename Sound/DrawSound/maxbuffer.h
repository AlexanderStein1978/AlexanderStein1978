#pragma once


#include "patternbuffer.h"

#include <QString>

class MaxBuffer
{
public:
    enum Observation {NothingObserved, BObserved, DObserved};

    constexpr static const double Invalid = -1.0;

    MaxBuffer(const double diameter);
    ~MaxBuffer();

    double newValue(const double f, const double a);
    double searchBroadest(const double f, const double a, const double pf);
    void getBroadest(const int index, double& start, double& end, double& max, QString& ratioForLabel);
    Observation analyzeNewValue(const double f, const double a);
    double getObsStart() const;
    double getObsEnd() const;
    double getMaxAmplitude() const;
    void clearObs();

    inline PatternBuffer::PatternObservation* popObservation()
    {
        return mPatternBuffer.popObservation();
    }

private:

    struct element
    {
        double f, a;
        element *prev, *next;
    };

    struct feature
    {
        double fStart = Invalid, fEnd, maxA = 0.0, ratio;
        feature *next = nullptr;
    };

    const double mMinWidth = 5e-4;
    const double mThreshold = 5e-5;

    void DestroyChain(element* first, element* last);
    void AddToStore(element* first, element* last);

    element *mFirst, *mLast, *mFStore, *mLStore;
    feature *mCurrent, *mFeatures;
    const double mDiameter;
    PatternBuffer mPatternBuffer;
    double mPrevTime, mSecondLastAmp = 0.0, mLastAmp = 0.0;
};
