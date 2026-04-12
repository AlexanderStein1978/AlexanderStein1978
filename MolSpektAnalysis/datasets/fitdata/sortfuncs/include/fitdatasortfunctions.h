//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SORTFUNCTIONS_H
#define SORTFUNCTIONS_H


class TableViewWindowCore;


bool sortByProg(const TableViewWindowCore *const core, const int n, const int m);
bool sortIvJFreqF(const TableViewWindowCore *const core, const int n, const int m);
bool SortvsIvJ(const TableViewWindowCore *const core, const int n, const int m);
bool sortIefJFreq(const TableViewWindowCore *const core, const int n, const int m);
bool sortIefJvFreq(const TableViewWindowCore *const core, const int n, const int m);
bool sortIefJFreqv(const TableViewWindowCore *const core, const int n, const int m);
bool sortbyDeviation(const TableViewWindowCore *const core, const int n, const int m);
bool sortbyDevR(const TableViewWindowCore *const core, const int n, const int m);
bool sortforTFGS(const TableViewWindowCore *const core, const int n, const int m);
bool sortByElState(const TableViewWindowCore *const core, const int n, const int m);
bool sortForExtractNewOrChanged(const TableViewWindowCore* const core, const int n, const int m);
bool sortByIsoColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortBy_vColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByJColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortBy_vsColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByJsColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortBySourceColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByProgressionColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByFileColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByEnergyColumn(const TableViewWindowCore* const core, const int n, const int m);
bool sortByUncertaintyColumn(const TableViewWindowCore* const core, const int n, const int m);



#endif
