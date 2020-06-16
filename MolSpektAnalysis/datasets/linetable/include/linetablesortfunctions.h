//
// C++ Interface: sortfunctions
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LINETABLE_SORTFUNCTIONS_H
#define LINETABLE_SORTFUNCTIONS_H


class QTableWidget;
class TermEnergy;


bool sortUtIvJ(const QTableWidget *const Tab, const int n, const int m);
bool sortIJvP(const QTableWidget *const Tab, const int n, const int m);
bool sortIvPJ(const QTableWidget *const Tab, const int n, const int m);
bool sortFPInt(const QTableWidget *const Tab, const int n, const int m);
bool sortIJvFreq(const QTableWidget *const Tab, const int n, const int m);
bool sortByvs(const QTableWidget *const Tab, const int n, const int m);
bool sortForSPN(const QTableWidget *const Tab, const int n, const int m);
bool sortBySpectrum(const QTableWidget *const Tab, const int n, const int m);
bool sortfRemDoubl(const QTableWidget *const Tab, const int n, const int m);
bool sortByFrequency(const QTableWidget *const Tab, const int n, const int m);
bool sortByProgression(const QTableWidget *const Tab, const int n, const int m);
bool isnSPG(const QTableWidget *const Tab, const int n, const int m);

bool isnSPG(TermEnergy& T1, TermEnergy &T2);

#endif
