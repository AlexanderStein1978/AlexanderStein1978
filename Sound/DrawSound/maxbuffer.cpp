#include "maxbuffer.h"


MaxBuffer::MaxBuffer(const double diameter) : mFirst(nullptr), mLast(nullptr), mFStore(nullptr), mLStore(nullptr), mDiameter(diameter)
{
}

MaxBuffer::~MaxBuffer()
{
    DestroyChain(mFirst, mLast);
    DestroyChain(mFStore, mLStore);
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
        first->prev = mLStore->prev;
        first->prev->next = first;
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
            while (nextLarger->a <= a) nextLarger = nextLarger->prev);
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
            mFStore->prev = mLast->prev;
            mLast = mFStore;
            mFStore = mFStore->next;
            mFStore->prev = nullptr;
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
