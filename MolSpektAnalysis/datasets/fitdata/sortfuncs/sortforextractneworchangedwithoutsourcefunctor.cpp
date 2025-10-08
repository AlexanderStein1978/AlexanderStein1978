//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "sortforextractneworchangedwithoutsourcefunctor.h"
#include "fitdata.h"

#include <QTableWidget>


SortForExtractNewOrChangedWithoutSourceFunctor::SortForExtractNewOrChangedWithoutSourceFunctor(const QTableWidget * const i_Tab, const QString * const i_SourceOffsetNames,
                                                                                               const double * const i_SourceOffset, const int i_NSourceOffset)
    : m_Tab(i_Tab), m_SourceOffsetNames(i_SourceOffsetNames), m_SourceOffset(i_SourceOffset), m_NSourceOffset(i_NSourceOffset)
{
}

bool SortForExtractNewOrChangedWithoutSourceFunctor::operator ()(const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const QString S1(m_Tab->item(n, fdcLineElState)->text()), S2(m_Tab->item(m, fdcLineElState)->text());
    if (S1 < S2) return true;
    if (S1 > S2) return false;
    const int Iso1 = m_Tab->item(n, fdcIso)->text().toInt(), Iso2 = m_Tab->item(m, fdcIso)->text().toInt();
    if (Iso1 < Iso2) return true;
    if (Iso1 > Iso2) return false;
    const int v1 = m_Tab->item(n, fdcv)->text().toInt(), v2 = m_Tab->item(m, fdcv)->text().toInt();
    if (v1 < v2) return true;
    if (v1 > v2) return false;
    const int J1 = m_Tab->item(n, fdcJ)->text().toInt(), J2 = m_Tab->item(m, fdcJ)->text().toInt();
    if (J1 < J2) return true;
    if (J1 > J2) return false;
    const bool ef1 = (m_Tab->item(n, fdcJs)->text().toInt() == J1), ef2 = (m_Tab->item(m, fdcJs)->text().toInt() == J2);
    if (ef1 && !ef2) return true;
    if (!ef1 && ef2) return false;
    const bool abs1 = m_Tab->item(n, fdcUncert)->text().contains(QRegExp("[234]01")), abs2 = m_Tab->item(m, fdcUncert)->text().contains(QRegExp("[234]01"));
    if (abs1 && !abs2) return true;
    if (!abs1 && abs2) return false;
    double DeltaE1 = 0.0, DeltaE2 = 0.0;
    if (m_NSourceOffset > 0)
    {
        QString SourceName1 = m_Tab->item(n, fdcSource)->text(), SourceName2 = m_Tab->item(m, fdcSource)->text();
        for (int i=0; i < m_NSourceOffset; ++i)
        {
            if (m_SourceOffsetNames[i] == SourceName1) DeltaE1 = m_SourceOffset[i];
            if (m_SourceOffsetNames[i] == SourceName2) DeltaE2 = m_SourceOffset[i];
        }
    }
    const double E1 = m_Tab->item(n, fdcEnergy)->text().toDouble() - DeltaE1, E2 = m_Tab->item(m, fdcEnergy)->text().toDouble() - DeltaE2;
    if (E1 < E2) return true;
    return false;
}