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
    delete mCore;
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

void TableViewWindow::getTabDimensions(int& NRows, int& NCols)
{
    NRows = mCore->rowCount();
    NCols = mCore->columnCount();
}

void TableViewWindow::setTabDimensions(int NRows, int)
{
    mCore->setRowCount(NRows);
}

void TableViewWindow::setRowData(int Row, QString* Data)
{
    QStringList L;
    int NC = mCore->columnCount();
    L.reserve(NC);
    table->blockSignals(true);
    mCore->blockSignals(true);
    for (int n=0; n < NC; ++n) L.push_back(Data[n]);
    mCore->setRow(L, Row);
    mCore->blockSignals(false);
    table->blockSignals(false);
    Changed();
}

bool TableViewWindow::checkAllConnections(int)
{
    bool ret = true, changed = false;
    Spektrum *Spectrum = new Spektrum(MW);
    QString SpektPath, NewPath;
    QStringList CheckedPaths, Replacements;
    int n, i, N = mCore->rowCount();
    table->blockSignals(true);
    mCore->blockSignals(true);
    for (n=0; n < N; ++n) if ((SpektPath = mCore->getSourceFile(n)).indexOf('(') == -1)
    {
        if ((i = CheckedPaths.indexOf(SpektPath)) == -1)
        {
            if (Spectrum->readData(SpektPath)) NewPath = Spectrum->getFileName();
            else
            {
                NewPath = SpektPath;
                ret = false;
            }
            CheckedPaths << SpektPath;
            Replacements << NewPath;
        }
        else NewPath = Replacements[i];
        if (NewPath != SpektPath)
        {
            mCore->setSourceFile(n, NewPath);
            changed = true;
        }
    }
    mCore->blockSignals(false);
    table->blockSignals(false);
    if (changed) Changed();
    delete Spectrum;
    return ret;
}
