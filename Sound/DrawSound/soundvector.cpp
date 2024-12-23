#include "soundvector.h"

#include <cstring>
#include <cmath>


SoundVector::SoundVector() : mData(nullptr), mSize(0)
{
}


SoundVector::SoundVector(const int numElements) : mData(new double[numElements]), mSize(numElements)
{
    memset(mData, 0, sizeof(double) * numElements);
}

SoundVector::~SoundVector()
{
    if (nullptr != mData) delete[] mData;
}

SoundVector::SoundVector(const SoundVector& temp) : mData(temp.mData), mSize(temp.mSize)
{
    const_cast<SoundVector*>(&temp)->mData = nullptr;
    const_cast<SoundVector*>(&temp)->mSize = 0;
}

double& SoundVector::operator[](const int index)
{
    return mData[index];
}

const double& SoundVector::operator[](const int index) const
{
    return mData[index];
}


int SoundVector::getSize() const
{
    return mSize;
}

void SoundVector::resize(const int newSize)
{
    if (nullptr == mData) mData = new double[mSize = newSize];
    else
    {
        if (newSize > mSize)
        {
            double *newData = new double[newSize];
            for (int n=0; n < mSize; ++n) newData[n] = mData[n];
            memset(mData + mSize, 0, sizeof(double) * (newSize - mSize));
            delete[] mData;
            mData = newData;
        }
        mSize = newSize;
    }
}

SoundVector SoundVector::sigmoid() const
{
    SoundVector result(mSize+1);
    result[0] = 1.0;
    for (int n=0; n < mSize; ++n) result[n+1] = 1.0 / (1.0 * exp(mData[n]));
    return result;
}

int SoundVector::getIndexOfMax() const
{
    int maxIndex = 0;
    for (int n=1; n < mSize; ++n) if (mData[n] > mData[maxIndex]) maxIndex = n;
    return maxIndex;
}
