//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SORTFUNCTIONS_H
#define SORTFUNCTIONS_H


class FitDataCore;


bool sortByProg(const FitDataCore *const Tab, const int n, const int m);
bool sortIvJFreqF(const FitDataCore *const Tab, const int n, const int m);
bool SortvsIvJ(const FitDataCore *const Tab, const int n, const int m);
bool sortIefJFreq(const FitDataCore *const Tab, const int n, const int m);
bool sortIefJvFreq(const FitDataCore *const Tab, const int n, const int m);
bool sortIefJFreqv(const FitDataCore *const Tab, const int n, const int m);
bool sortbyDeviation(const FitDataCore *const Tab, const int n, const int m);
bool sortbyDevR(const FitDataCore *const Tab, const int n, const int m);
bool sortforTFGS(const FitDataCore *const Tab, const int n, const int m);
bool sortByElState(const FitDataCore *const Tab, const int n, const int m);
bool sortForExtractNewOrChanged(const FitDataCore* const Tab, const int n, const int m);

#endif
