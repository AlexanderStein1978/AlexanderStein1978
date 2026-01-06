//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2026
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LINETABLE_SORTFUNCTIONS_H
#define LINETABLE_SORTFUNCTIONS_H


class TableViewWindowCore;
class TermEnergy;


bool sortUtIvJ(const TableViewWindowCore *const core, const int n, const int m);
bool sortIJvP(const TableViewWindowCore *const core, const int n, const int m);
bool sortIvPJ(const TableViewWindowCore *const core, const int n, const int m);
bool sortFPInt(const TableViewWindowCore *const core, const int n, const int m);
bool sortIJvFreq(const TableViewWindowCore *const core, const int n, const int m);
bool sortByvs(const TableViewWindowCore *const core, const int n, const int m);
bool sortForSPN(const TableViewWindowCore *const core, const int n, const int m);
bool sortBySpectrum(const TableViewWindowCore *const core, const int n, const int m);
bool sortfRemDoubl(const TableViewWindowCore *const core, const int n, const int m);
bool sortByFrequency(const TableViewWindowCore *const core, const int n, const int m);
bool sortByProgression(const TableViewWindowCore *const core, const int n, const int m);
bool isnSPG(const TableViewWindowCore *const core, const int n, const int m);

bool isnSPG(TermEnergy& T1, TermEnergy &T2);

#endif
