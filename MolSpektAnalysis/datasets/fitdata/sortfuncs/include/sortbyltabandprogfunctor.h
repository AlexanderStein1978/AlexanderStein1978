//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SORTBYLTABANDPROGFUNCTOR_H
#define SORTBYLTABANDPROGFUNCTOR_H


class QTableWidget;

class LineTable;
class MainWindow;
class Molecule;


class SortByLTabAndProgFunctor
{
public:
    SortByLTabAndProgFunctor(QTableWidget* Tab, LineTable** Sources, int NSources, MainWindow *MW, Molecule *Mol) : m_Tab(Tab), m_Sources(Sources), m_NSources(NSources), m_MW(MW), m_Mol(Mol)
    {}

    bool operator()(int n, int m);

private:
    QTableWidget* m_Tab;
    LineTable** m_Sources;
    int m_NSources;
    MainWindow* m_MW;
    Molecule* m_Mol;
};

#endif
