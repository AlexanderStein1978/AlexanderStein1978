#include "datasortfunctor.h"
#include "utils.h"


DataSortFunctor::DataSortFunctor(double ** data, int numRows, int numColumns) : mNumData(numRows * numColumns), mHelper(CreateInt(mNumData, 2)), mData(data)
{
    int n=-1;
    for (int r=0; r < numRows; ++r) for (int c=0; c < numColumns; ++c)
    {
        mHelper[++n][0] = r;
        mHelper[n][1] = c;
    }
}

bool DataSortFunctor::operator()(const int index1, const int index2) const
{
    if (index1 == -1) return false;
    if (index2 == -1) return true;
    return mData[mHelper[index1][0]][mHelper[index1][1]] < mData[mHelper[index2][0]][mHelper[index2][1]];
}

int ** DataSortFunctor::getResult(int* sortOrder)
{
    int** result = new int*[mNumData];
    for (int n=0; n < mNumData; ++n) result[sortOrder[n]] = mHelper[n];
    delete[] sortOrder;
    delete[] mHelper;
    mHelper = nullptr;
    return result;
}
