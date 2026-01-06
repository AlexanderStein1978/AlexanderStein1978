//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2026
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LINETABLE_SORTFUNCTIONS_H
#define LINETABLE_SORTFUNCTIONS_H


class LineTableCore;
class TermEnergy;


bool sortUtIvJ(const LineTableCore *const core, const int n, const int m);
bool sortIJvP(const LineTableCore *const core, const int n, const int m);
bool sortIvPJ(const LineTableCore *const core, const int n, const int m);
bool sortFPInt(const LineTableCore *const core, const int n, const int m);
bool sortIJvFreq(const LineTableCore *const core, const int n, const int m);
bool sortByvs(const LineTableCore *const core, const int n, const int m);
bool sortForSPN(const LineTableCore *const core, const int n, const int m);
bool sortBySpectrum(const LineTableCore *const core, const int n, const int m);
bool sortfRemDoubl(const LineTableCore *const core, const int n, const int m);
bool sortByFrequency(const LineTableCore *const core, const int n, const int m);
bool sortByProgression(const LineTableCore *const core, const int n, const int m);
bool isnSPG(const LineTableCore *const core, const int n, const int m);

bool isnSPG(TermEnergy& T1, TermEnergy &T2);

#endif
