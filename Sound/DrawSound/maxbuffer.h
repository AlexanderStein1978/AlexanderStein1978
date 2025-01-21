#pragma once


class MaxBuffer
{
public:
    MaxBuffer(const double diameter);
    ~MaxBuffer();

    double newValue(const double f, const double a);

private:

    struct element
    {
        double f, a;
        element *prev, *next;
    };

    void DestroyChain(element* first, element* last);
    void AddToStore(element* first, element* last);

    element* mFirst, *mLast, *mFStore, *mLStore;
    const double mDiameter;
};
