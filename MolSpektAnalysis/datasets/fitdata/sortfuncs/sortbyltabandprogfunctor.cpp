//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "sortbyltabandprogfunctor.h"
#include "MainWindow.h"
#include "fitdata.h"

#include <QTableWidget>


bool SortByLTabAndProgFunctor::operator()(int n, int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    LineTable* LineTabn = (n < m_NSources ? m_Sources[n] : m_MW->getLineTable(m_Tab->item(n, fdcSource)->text(), m_Mol));
    LineTable* LineTabm = (m < m_NSources ? m_Sources[m] : m_MW->getLineTable(m_Tab->item(m, fdcSource)->text(), m_Mol));
    if (LineTabn < LineTabm) return true;
    if (LineTabn > LineTabm) return false;
    if (m_Tab->item(n, fdcProg)->text().toInt() < m_Tab->item(m, fdcProg)->text().toInt()) return true;
    return false;
}