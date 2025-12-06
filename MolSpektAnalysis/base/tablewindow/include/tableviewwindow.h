//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_H
#define TABLE_VIEW_WINDOW_H


#include "tablewindow.h"


class TableViewWindowCore;


class TableViewWindow : public TableWindow
{
public:
    TableViewWindow(TableViewWindowCore *const mCore, Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    ~TableViewWindow();

    int getMaxJ();
    int getMaxv();
    bool isDataAvailable();

protected:
    TableViewWindowCore *const mCore;
};

#endif
