//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once


#include "roundbuffer.h"


class PatternBuffer : public RoundBuffer
{
public:
    enum ObservationStatus{NoPatternObserved, MaybePatternObserved, PatternObserved};

    struct Pattern
    {
        constexpr static const int  MaxLength = 5;

        Pattern();

        RoundBuffer buffer;
        Pattern* next;
        int length, index;
        double timeLength;
    };

    struct PatternObservation
    {
        int patternIndex;
        double start, end, maxA = 0.0;
        PatternObservation* next = nullptr;
    };

    PatternBuffer(const int size);
    virtual ~PatternBuffer();

    void setNextElement(const double frequency, const double amplitude) override;
    PatternObservation* popObservation();

private:
    const double AmplitudeTol = 0.1, TimeTol = 0.1;

    ObservationStatus checkPattern(Pattern* pattern, int startIndex);
    void createPattern();
    void createObservation();

    Pattern* mPatterns, *mLastPattern;
    Pattern* mCurrentPattern, mPossiblePattern, *mMaybeReobservedPattern = nullptr;
    PatternObservation* mObservations = nullptr;
    double mPossibleObsStart, mPossibleObsEnd, mMax = 0.0;
    ObservationStatus mObsStatus = NoPatternObserved;
};
