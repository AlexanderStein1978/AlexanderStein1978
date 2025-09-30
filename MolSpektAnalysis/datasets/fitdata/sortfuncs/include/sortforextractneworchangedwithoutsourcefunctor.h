//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SORTFOREXTRACTNEWORCHANGEDWITHOUTSOURCEFUNCTOR_H
#define SORTFOREXTRACTNEWORCHANGEDWITHOUTSOURCEFUNCTOR_H


class QTableWidget;
class QString;


class SortForExtractNewOrChangedWithoutSourceFunctor
{
public:
    SortForExtractNewOrChangedWithoutSourceFunctor(const QTableWidget * const i_Tab, const QString* const i_SourceOffsetNames, const double* const i_SourceOffset,
                                                   const int i_NSourceOffset);
    bool operator()(const int n, const int m);

private:
    const QTableWidget* const m_Tab;
    const QString* const m_SourceOffsetNames;
    const double* const m_SourceOffset;
    const int m_NSourceOffset;
};

#endif