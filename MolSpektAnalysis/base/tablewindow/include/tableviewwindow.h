//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_H
#define TABLE_VIEW_WINDOW_H


#include "tablewindow.h"


class TableViewWindow : public TableWindow
{
public:
    TableViewWindow(Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    ~TableViewWindow();
};

#endif
