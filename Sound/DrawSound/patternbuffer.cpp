#include "patternbuffer.h"

#include <cmath>


PatternBuffer::PatternBuffer(const int size) : RoundBuffer(size), mPatterns(nullptr), mLastPattern(nullptr), mCurrentPattern(nullptr), mObservations(nullptr), mPossibleObsStart(-1.0),
    mPossibleObsEnd(-1.0)
{
}

PatternBuffer::~PatternBuffer() noexcept
{
    while (nullptr != mPatterns)
    {
        Pattern* d = mPatterns;
        mPatterns = d->next;
        delete d;
    }
    while (nullptr != mObservations) delete popObservation();
}

void PatternBuffer::setNextElement(const double frequency, const double amplitude)
{
    mPossiblePattern.buffer.setNextElement(element(0).frequency - frequency, amplitude / element(0).amplitude);
    RoundBuffer::setNextElement(frequency, amplitude);
    int NElements = getNumElements();
    if (NElements < 3) return;
    if (nullptr != mCurrentPattern)
    {
        if (amplitude > mMax) mMax = amplitude;
        if ((frequency - mPossibleObsEnd) * (1.0 + TimeTol) < mCurrentPattern->timeLength) return;
        mObsStatus = checkPattern(mCurrentPattern, 0);
        if (mObsStatus == MaybePatternObserved) mObsStatus = PatternObserved;
        if (PatternObserved == mObsStatus) mPossibleObsEnd = frequency;
    }
    else if (mObsStatus == MaybePatternObserved)
    {
        if (nullptr == mMaybeReobservedPattern) mObsStatus = checkPattern(&mPossiblePattern, 0);
        else mObsStatus = checkPattern(mMaybeReobservedPattern, 0);
        if (MaybePatternObserved == mObsStatus)
        {
            mObsStatus = PatternObserved;
            if (nullptr == mMaybeReobservedPattern) createPattern();
            else
            {
                createObservation();
                mCurrentPattern = mMaybeReobservedPattern;
                mMaybeReobservedPattern = nullptr;
            }
            mPossibleObsStart = element(2).frequency;
            mPossibleObsEnd = frequency;
            mMax = (element(2).amplitude > element(1).amplitude ? element(2).amplitude : element(1).amplitude);
            if (amplitude > mMax) mMax = amplitude;
        }
    }
    if (mObsStatus == NoPatternObserved) for (Pattern* p = mPatterns; nullptr != p; p = p->next) if (NElements >= p->length)
    {
        mObsStatus = checkPattern(p, 0);
        if (mObsStatus == PatternObserved || mObsStatus == MaybePatternObserved)
        {
            if (mObsStatus == PatternObserved)
            {
                createObservation();
                mCurrentPattern = p;
                mPossibleObsStart = element(p->length).frequency;
                mPossibleObsEnd = frequency;
                mMax = amplitude;
                for (int i=1; i <= p->length; ++i) if (element(i).amplitude > mMax) mMax = element(i).amplitude;
            }
            else mMaybeReobservedPattern = p;
            break;
        }
    }
    if (mObsStatus == NoPatternObserved) for (int n=2; n <= mPossiblePattern.MaxLength && n < NElements; ++n)
    {
        mPossiblePattern.length = n;
        mObsStatus = checkPattern(&mPossiblePattern, n);
        if (mObsStatus == PatternObserved)
        {
            createObservation();
            createPattern();
            mPossibleObsEnd = frequency;
            mMax = amplitude;
            for (int i=1; i <= n; ++i) if (element(i).amplitude > mMax) mMax = element(i).amplitude;
            break;
        }
    }
    if (mObsStatus == NoPatternObserved)
    {
        mPossiblePattern.length = 1;
        mObsStatus = checkPattern(&mPossiblePattern, 1);
        if (mObsStatus == MaybePatternObserved)
        {
            mObsStatus = checkPattern(&mPossiblePattern, 2);
            if (mObsStatus = NoPatternObserved) mObsStatus = MaybePatternObserved;
            else
            {
                mObsStatus = PatternObserved;
                createObservation();
                createPattern();
                mPossibleObsStart = element(2).frequency;
                mPossibleObsEnd = frequency;
                mMax = (element(2).amplitude > element(1).amplitude ? element(2).amplitude : element(1).amplitude);
                if (amplitude > mMax) mMax = amplitude;
            }
        }
    }
}

PatternBuffer::ObservationStatus PatternBuffer::checkPattern(Pattern* pattern, int startIndex)
{
    for (int n=0; n < pattern->length; ++n)
    {
        double expT = element(startIndex + n).frequency + pattern->buffer.element(n).frequency;
        double expA = element(startIndex + n).amplitude * pattern->buffer.element(n).amplitude;
        if (abs(element(startIndex + n + 1).frequency - expT) > TimeTol * (element(startIndex + n + 1).frequency - element(startIndex + n).frequency)
            || abs(element(startIndex + n).amplitude - expA) > AmplitudeTol * element(startIndex + n).amplitude)
        {
            return NoPatternObserved;
        }
    }
    return (pattern->length > 1 ? PatternObserved : MaybePatternObserved);
}

void PatternBuffer::createPattern()
{
    createObservation();
    mCurrentPattern = new Pattern(mPossiblePattern);
    mCurrentPattern->index = mLastPattern->index + 1;
    mLastPattern->next = mCurrentPattern;
    mLastPattern = mCurrentPattern;
}

void PatternBuffer::createObservation()
{
    if (nullptr == mCurrentPattern) return;
    if (mObservations->patternIndex == mCurrentPattern->index)
    {
        mObservations->end = mPossibleObsEnd;
        return;
    }
    PatternObservation* newObservation = new PatternObservation;
    newObservation->end = mPossibleObsEnd;
    newObservation->maxA = mMax;
    newObservation->patternIndex = mCurrentPattern->index;
    newObservation->start = mPossibleObsStart;
    newObservation->next = mObservations;
    mObservations = newObservation;
}

PatternBuffer::PatternObservation* PatternBuffer::popObservation()
{
    PatternObservation* r = mObservations;
    mObservations = r->next;
    return r;
}

PatternBuffer::Pattern::Pattern() : buffer(MaxLength), next(nullptr), index(0), length(0), timeLength(0.0)
{
}

