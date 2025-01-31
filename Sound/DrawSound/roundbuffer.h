#pragma once


class RoundBuffer
{
public:
    enum ObservationStatus{NoPatternObserved, MaybePatternObserved, PatternObserved};
    enum PatternStatus{NewPattern, ReobservedPattern};

    struct BufferElement
    {
        double frequency, amplitude;
    };

    struct Pattern
    {
        constexpr static const int  MaxLength = 5;

        BufferElement elements[MaxLength];
        Pattern* next;
        int length = 0;
    };

    struct PatternObservation
    {
        int patternIndex;
        double start, end, maxA = 0.0;
        PatternObservation* next = nullptr;
    };

    RoundBuffer(const int size);
    ~RoundBuffer();

    void setNextElement(const double frequency, const double amplitude);
    void setElementToAnalyse(const double frequency, const double amplitude);
    const BufferElement& element(const int index) const;
    PatternObservation* popObservation();

private:
    BufferElement* m_elements;
    Pattern* mPatterns, *mLastPattern;
    Pattern* mCurrentPattern;
    PatternObservation* observations;
    int m_size, m_currentIndex, mNumElements;
    double mPossibleObsStart, mPossibleObsEnd;
    ObservationStatus mObsStatus;
    PatternStatus mPatternStatus;
};
