//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2026
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLEVEWWINDOWCORESORTFUNCTOR_H
#define TABLEVEWWINDOWCORESORTFUNCTOR_H


class TableViewWindowCore;


class TableViewWindowCoreSortFunctor
{
public:
    TableViewWindowCoreSortFunctor(TableViewWindowCore* core, bool sortFuncs(const TableViewWindowCore* const, const int, const int))
        : mCore(core), m_sortFuncs(*sortFuncs)
    {
    }

    inline bool operator()(int i, int j)
    {
        return m_sortFuncs(mCore, i, j);
    }

private:
    TableViewWindowCore* mCore;
    bool (&m_sortFuncs)(const TableViewWindowCore* const, const int, const int);
};

#endif
