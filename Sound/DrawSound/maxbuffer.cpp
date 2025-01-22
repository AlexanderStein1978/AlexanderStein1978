#include "maxbuffer.h"

#include <cmath>


MaxBuffer::MaxBuffer(const double diameter) : mFirst(nullptr), mLast(nullptr), mFStore(nullptr), mLStore(nullptr), mCurrent(new feature), mFeatures(nullptr), mDiameter(diameter)
{
}

MaxBuffer::~MaxBuffer()
{
    DestroyChain(mFirst, mLast);
    DestroyChain(mFStore, mLStore);
    clearObs();
    delete mCurrent;
}

void MaxBuffer::DestroyChain(element* first, element* last)
{
    if (nullptr == first) return;
    for (element* c = first->next; nullptr != c; c = c->next) delete c->prev;
    delete last;
}

void MaxBuffer::AddToStore(element* first, element* last)
{
    if (nullptr == first) return;
    last->next = nullptr;
    if (nullptr == mFStore)
    {
        mFStore = first;
        first->prev = nullptr;
    }
    else
    {
        first->prev = mLStore;
        mLStore->next = first;
    }
    mLStore = last;
}

double MaxBuffer::newValue(const double f, const double a)
{
    double prevMax = 0.0;;
    if (nullptr == mFirst)
    {
        mFirst = mLast = new element;
        mFirst->prev = nullptr;
    }
    else
    {
        prevMax = mFirst->a;
        if (f - mFirst->f > mDiameter)
        {
            if (mFirst == mLast) mFirst->a = 0.0;
            else
            {
                element* newFirst = mFirst->next;
                AddToStore(mFirst, mFirst);
                mFirst = newFirst;
                mFirst->prev = nullptr;
            }
        }
        if (a >= mFirst->a)
        {
            AddToStore(mFirst->next, mLast);
            mLast = mFirst;
        }
        else if (a >= mLast->a)
        {
            element *nextLarger = mLast;
            while (nextLarger->a <= a) nextLarger = nextLarger->prev;
            if (nextLarger->next != mLast)
            {
                element* current = nextLarger->next;
                AddToStore(current->next, mLast);
                mLast = current;
            }
        }
        else if (nullptr != mFStore)
        {
            mLast->next = mFStore;
            mFStore->prev = mLast;
            mLast = mFStore;
            mFStore = mFStore->next;
            if (nullptr != mFStore) mFStore->prev = nullptr;
            else mLStore = nullptr;
        }
        else
        {
            mLast->next = new element;
            mLast->next->prev = mLast;
            mLast = mLast->next;
        }
    }
    mLast->next = nullptr;
    mLast->a = a;
    mLast->f = f;
    return mFirst->a - prevMax;
}

MaxBuffer::Observation MaxBuffer::analyzeNewValue(const double f, const double a)
{
    double newVal = newValue(f, a);
    if (newVal > mThreshold)
    {
        double Abs = fabs(a);
        if (Abs > mCurrent->maxA) mCurrent->maxA = Abs;
        if (mCurrent->fStart == Invalid) mCurrent->fStart = f;
        mCurrent->fEnd = f;
    }
    else if (Invalid != mCurrent->fStart && f - mCurrent->fEnd > mMinWidth)
    {
        if (mCurrent->fEnd - mCurrent->fStart < mMinWidth)
        {
            mCurrent->maxA = 0.0;
            mCurrent->fStart = Invalid;
            if (nullptr != mFeatures)
            {
                if (nullptr == mFeatures->next) return DObserved;
                else return BObserved;
            }
        }
        else
        {
            mCurrent->next = mFeatures;
            mFeatures = mCurrent;
            mCurrent = new feature;
        }
    }
    return NothingObserved;
}

void MaxBuffer::clearObs()
{
    if (nullptr != mFeatures)
    {
        feature *c, *n;
        for (c = mFeatures; c != nullptr; c=n)
        {
            n = c->next;
            delete c;
        }
        mFeatures = nullptr;
    }
}

double MaxBuffer::getMaxAmplitude() const
{
    double result = Invalid;
    for (feature* c = mFeatures; c != nullptr; c = c->next) if (c->maxA > result) result = c->maxA;
    return result;
}

double MaxBuffer::getObsEnd() const
{
    if (mFeatures != nullptr) return mFeatures->fEnd;
    return Invalid;
}

double MaxBuffer::getObsStart() const
{
    if (nullptr == mFeatures) return Invalid;
    feature* f = mFeatures;
    while (nullptr != f->next) f = f->next;
    return f->fStart;
}
