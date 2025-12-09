//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tableviewwindow.h"
#include "tableviewwindowcore.h"
#include "Spektrum.h"
#include "utils.h"

#include <QTextStream>


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

void TableViewWindow::copyRows(int& numRows, int& numColums, int *& Rows, QString **& Data)
{
	if (Data != 0) Destroy(Data, numRows);
	table->getSelectedRows(Rows, numRows);
    numColums = mCore->columnCount();
	Data = CreateQString(numRows, numColums);
	for (int r=0; r <  numRows; ++r) BaseDataToQStringArray(*mCore->getRow(Rows[r]), Data[r]);
}

void TableViewWindow::copyRows(int& numRows, int& numColums, QString **& Data)
{
	int *Rows;
	copyRows(numRows, numColums, Rows, Data);
    delete[] Rows;
    Changed();
}

void TableViewWindow::cutRows(int& numRows, int& numColums, QString **& Data)
{
    int *Rows;
	copyRows(numRows, numColums, Rows, Data);
    table->blockSignals(true);
    mCore->blockSignals(true);
    mCore->deleteRows(Rows, numRows);
    mCore->blockSignals(false);
    table->blockSignals(false);
    delete[] Rows;
    Changed();
}

void TableViewWindow::insertRows(int numRows, int numColumns, QString ** Data)
{
    table->blockSignals(true);
    mCore->blockSignals(true);
    for (int r=0; r < numRows; ++r)
    {
        QStringList L;
        for (int c=0; c < numColumns; ++c) L << Data[r][c];
        mCore->addRow(L);
    }
    mCore->blockSignals(false);
    table->blockSignals(false);
    Changed();
}

void TableViewWindow::MarkLines(int* rN, int N)
{
    int n=0, r1, MC = mCore->columnCount() - 1;
	QItemSelectionModel* model = new QItemSelectionModel;
    QModelIndex topIndex = mCore->getIndex(rN[0], 0);
    table->blockSignals(true);
	table->scrollTo(topIndex, QAbstractItemView::PositionAtTop);
	for (r1 = 1; r1 != 0; ) for (r1 = 0, n=1; n<N; n++) if (rN[n] < rN[n-1])
	{
		r1 = rN[n];
		rN[n] = rN[n-1];
		rN[n-1] = r1;
	}
	for (n=0; n<N; )
	{
		for (r1 = rN[n++]; n<N && rN[n] == rN[n-1] + 1; ++n) ;
        QModelIndex topLeft = mCore->getIndex(r1, 0), bottomRight = mCore->getIndex(rN[n-1], MC);
        QItemSelection newSelection(topLeft, bottomRight);
        model->select(newSelection, QItemSelectionModel::Select);
	}
	table->setSelectionModel(model);
    table->blockSignals(false);
	if (!isVisible()) show();
	activateWindow();
	table->setFocus();
}

void TableViewWindow::exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns)
{
	int n, r, c, NR = mCore->rowCount(), NC = mCore->columnCount();
	QFile F(FileName);
	F.open(QIODevice::WriteOnly);
	QTextStream S(&F);
	if (selectedCells)
	{
        QModelIndexList L = table->getSelectedIndexes();
		bool *verticalCells = new bool[NR], *horizontalCells = new bool[NC], **cells = CreateBool(NR, NC);
        int top, bottom, left, right;
        for (n=0; n < NR; ++n) verticalCells[n] = false;
        for (n=0; n < NC; ++n) horizontalCells[n] = false;
        for (r=0; r < NR; ++r) for (c=0; c < NC; ++c) cells[r][c] = false;
        for (auto it = L.begin(); it != L.end(); ++it)
        {
            r = it->row();
            c = it->column();
            verticalCells[r] = true;
            horizontalCells[c] = true;
            cells[r][c] = true;
        }
        for (top=0; top < NR; ++top) if (verticalCells[top]) break;
        for (bottom = NR - 1; bottom >= 0; --bottom) if (verticalCells[bottom]) break;
        for (left=0; left < NC; ++left) if (horizontalCells[left]) break;
        for (right = NC - 1; right >= 0; --right) if (horizontalCells[right]) break;
		if (exchangeRowsColumns)
		{
			for (r = top; r <= bottom; ++r) if (verticalCells[r]) S << '\t' << mCore->headerData(r, Qt::Vertical).toByteArray();
			S << '\n';
			for (c = left; c <= right; ++c)
			{
				S << mCore->headerData(c, Qt::Horizontal).toByteArray();
				for (r = top; r <= bottom; ++r)
                {
                    S  << '\t';
					if (cells[r][c]) S << mCore->cellToString(r, c);
                }
                S << (r < bottom ? '\t' : '\n');
			}
		}
		else for (n=0; n < L.count(); ++n)
		{
			for (c = left; c <= right; ++c) S << '\t' << mCore->headerData(c, Qt::Horizontal).toByteArray();
			S << '\n';
			for (r = top; r <= bottom; ++r)
			{
				S << mCore->headerData(r, Qt::Vertical).toByteArray();
				for (c = left; c <= right; ++c)
                {
                    S << '\t';
                    if (cells[r][c]) S << mCore->cellToString(r, c);
                }
                S  << (c < right ? '\t' : '\n');
			}
		}
		delete[] horizontalCells;
        delete[] verticalCells;
        Destroy(cells, NR);
	}
	else
	{
		if (exchangeRowsColumns)
		{
			for (r=0; r < NR; ++r)	S << '\t' << mCore->headerData(r, Qt::Vertical).toByteArray();
			S << '\n';
			for (c=0; c < NC; ++c)
			{
				S << mCore->headerData(c, Qt::Horizontal).toByteArray() << '\t';
				for (r=0; r < NR; ++r) S << mCore->cellToString(r, c)  << (r < NR - 1 ? '\t' : '\n');
			}
		}
		else
		{
			for (c=0; c < NC; c++) S << '\t' << mCore->headerData(c, Qt::Horizontal).toByteArray();
			S << '\n';
			for (r=0; r < NR; r++)
			{
				S << mCore->headerData(r, Qt::Vertical).toByteArray() << '\t';
				for (c=0; c < NC; c++) S  << mCore->cellToString(r, c) << (c < NC - 1 ? '\t' : '\n');
			}
		}
	}
}

QString ** TableViewWindow::getData(int& NRows, int& NCols)
{
    int r, NR = mCore->rowCount(), NC = mCore->columnCount();
	if (NRows <= 0 || NRows > NR) NRows = NR;
	if (NCols <= 0 || NCols > NC) NCols = NC;
	QString **Data = CreateQString(NRows, NCols);
	for (r=0; r < NRows; ++r) BaseDataToQStringArray(*mCore->getRow(r), Data[r]);
	return Data;
}

void TableViewWindow::ContentChanged(const QModelIndex& topLeft, const QModelIndex& bottomRight, const QVector<int>&)
{
    int startRow = topLeft.row(), endRow = bottomRight.row(), r, N = mCore->columnCount();
	QString *IT = new QString[N];
    for (r = startRow; r <= endRow; ++r)
    {
        BaseDataToQStringArray(*mCore->getRow(r), IT);
        emit TabRowChanged(IT, N);
    }
	delete[] IT;
}

void TableViewWindow::setData(QString ** Data, int NRows, int NCols)
{
    int r, c, NC = mCore->columnCount();
    QStringList L;
	table->blockSignals(true);
    mCore->blockSignals(true);
	mCore->setRowCount(NRows);
	for (r=0; r < NRows; r++)
    {
        for (c=0; c < NC; c++)
        {
            if (c < NCols) L << Data[r][c];
            else L << "";
        }
        mCore->setRow(mCore->convertToBaseData(L), r);
        L.clear();
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}
