#pragma once


class MaxBuffer
{
public:
    enum Observation {NothingObserved, BObserved, DObserved};

    const double Invalid = -1.0;

    MaxBuffer(const double diameter);
    ~MaxBuffer();

    double newValue(const double f, const double a);
    Observation analyzeNewValue(const double f, const double a);
    double getObsStart() const;
    double getObsEnd() const;
    double getMaxAmplitude() const;
    void clearObs();

private:

    struct element
    {
        double f, a;
        element *prev, *next;
    };

    struct feature
    {
        double fStart = -1.0, fEnd, maxA = 0.0;
        feature *next = nullptr;
    };

    const double mMinWidth = 5e-4;
    const double mThreshold = 5e-5;

    void DestroyChain(element* first, element* last);
    void AddToStore(element* first, element* last);

    element *mFirst, *mLast, *mFStore, *mLStore;
    feature *mCurrent, *mFeatures;
    const double mDiameter;
};
