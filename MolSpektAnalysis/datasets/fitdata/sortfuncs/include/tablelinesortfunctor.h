//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TABLELINESORTFUNCTOR_H
#define TABLELINESORTFUNCTOR_H


class TableLine;


class TableLineSortFunctor
{
public:
    inline TableLineSortFunctor(TableLine** tableLinesToSort, int *NTableLinesForState)
        : TLToSort(tableLinesToSort), NTLS(NTableLinesForState)
    {
    }

    bool operator()(int i, int j);

private:
    TableLine** TLToSort;
    int *NTLS;
};

#endif
