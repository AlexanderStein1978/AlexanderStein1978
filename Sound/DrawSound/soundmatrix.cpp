#include "soundmatrix.h"
#include "soundvector.h"


SoundMatrix::SoundMatrix(const int numColumns, const int numRows) : mData(new SoundVector[numColumns]), mColumnCount(numColumns), mRowCount(numRows)
{
    for (int n=0; n < numColumns; ++n) mData[n].resize(numRows);
}

SoundMatrix::~SoundMatrix()
{
    delete[] mData;
}

SoundVector& SoundMatrix::operator[](const int column)
{
    return mData[column];
}

SoundVector SoundMatrix::operator*(const SoundVector& right)
{
    if (right.getSize() != mColumnCount) return SoundVector();
    SoundVector result(mRowCount);
    for (int c=0; c < mColumnCount; ++c) for (int r=0; r < mRowCount; ++r) result[r] += mData[c][r] * right[c];
    return result;
}
