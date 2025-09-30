//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

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
