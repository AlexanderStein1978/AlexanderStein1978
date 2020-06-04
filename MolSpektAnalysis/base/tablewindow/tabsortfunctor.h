//
// C++ Interface: TabSortFunctor
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABSORTFUNCTOR_H
#define TABSORTFUNCTOR_H


class TabSortFunctor
{
public:
    TabSortFunctor(QTableWidget* Tab, bool sortFuncs(const QTableWidget* const, const int, const int))
        : m_Tab(Tab), m_sortFuncs(*sortFuncs)
    {
    }

    inline bool operator()(int i, int j)
    {
        return m_sortFuncs(m_Tab, i, j);
    }

private:
    QTableWidget* m_Tab;
    bool (&m_sortFuncs)(const QTableWidget* const, const int, const int);
};

#endif
