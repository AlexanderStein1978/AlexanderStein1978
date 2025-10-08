//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SORTFUNCTIONS_H
#define SORTFUNCTIONS_H


class QTableWidget;


bool sortByProg(const QTableWidget *const Tab, const int n, const int m);
bool sortIvJFreqF(const QTableWidget *const Tab, const int n, const int m);
bool SortvsIvJ(const QTableWidget *const Tab, const int n, const int m);
bool sortIefJFreq(const QTableWidget *const Tab, const int n, const int m);
bool sortIefJvFreq(const QTableWidget *const Tab, const int n, const int m);
bool sortIefJFreqv(const QTableWidget *const Tab, const int n, const int m);
bool sortbyDeviation(const QTableWidget *const Tab, const int n, const int m);
bool sortbyDevR(const QTableWidget *const Tab, const int n, const int m);
bool sortforTFGS(const QTableWidget *const Tab, const int n, const int m);
bool sortByElState(const QTableWidget *const Tab, const int n, const int m);
bool sortForExtractNewOrChanged(const QTableWidget* const Tab, const int n, const int m);

#endif
