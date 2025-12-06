//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tableviewwindow.h"
#include "tableviewwindowcore.h"


TableViewWindow::TableViewWindow(TableViewWindowCore *const core, Type typ, MainWindow *MW, Molecule *M) : TableWindow(typ, MW, M), mCore(core)
{
}

TableViewWindow::~TableViewWindow()
{
}


int TableViewWindow::getMaxJ()
{
    return mCore->getMaxJ();
}

int TableViewWindow::getMaxv()
{
    return mCore->getMaxv();
}

bool TableViewWindow::isDataAvailable()
{
    if (mCore->rowCount() > 0) return true;
    return false;
}
