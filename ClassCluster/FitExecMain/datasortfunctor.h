#pragma once


class DataSortFunctor
{
public:
    DataSortFunctor(double **data, int numRows, int numColumns);

    bool operator()(const int index1, const int index2) const;

    int** getResult(int* sortOrder);

private:
    int mNumData;
    int** mHelper;
    double **mData;
};
