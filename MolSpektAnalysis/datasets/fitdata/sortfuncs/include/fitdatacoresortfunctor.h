//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FITDATACORESORTFUNCTOR_H
#define FITDATACORESORTFUNCTOR_H


class FitDataCore;


class FitDataCoreSortFunctor
{
public:
    FitDataCoreSortFunctor(FitDataCore* Tab, bool sortFuncs(const FitDataCore* const, const int, const int))
        : m_Tab(Tab), m_sortFuncs(*sortFuncs)
    {
    }

    inline bool operator()(int i, int j)
    {
        return m_sortFuncs(m_Tab, i, j);
    }

private:
    FitDataCore* m_Tab;
    bool (&m_sortFuncs)(const FitDataCore* const, const int, const int);
};

#endif
