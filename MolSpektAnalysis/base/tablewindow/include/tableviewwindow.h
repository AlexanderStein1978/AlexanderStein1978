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

struct BaseData;


class TableViewWindow : public TableWindow
{
    Q_OBJECT

public:
    TableViewWindow(TableViewWindowCore *const mCore, Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    virtual ~TableViewWindow();

    int getMaxJ();
    int getMaxv();
    bool isDataAvailable();
    void setTabDimensions(int NRows, int NCols) override;
	void getTabDimensions(int &NRows, int &NCols) override;
    void setRowData(int Row, QString *Data) override;
    void copyRows(int &numRows, int &numColums, QString **&Data) override;
	void cutRows(int &numRows, int &numColums, QString **&Data) override;
	void insertRows(int numRows, int numColumns, QString **Data) override;
	void MarkLines(int *rN, int N) override;
    void exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns) override;
    QString **getData(int &NRows, int &NCols) override;
    void setData(QString **Data, int NRows, int NCols) override;
    void scrollTo(const int row);

protected slots:
    void ContentChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight, const QVector<int> &roles = QVector<int>());
    void DeselectEverything();

	void HeaderItemClicked(const int index)
	{
		lastClickedHeaderIndex = index;
	}

protected:
    void resizeHelper(QRect& G) override;
    void writeData(QTextStream& S) override;
    void copyRows(int &numRows, int &numColums, int *&Rows, QString **&Data);
    virtual void BaseDataToQStringArray(const BaseData&, QString *const) const = 0;
    bool checkAllConnections(int FileColumn) override;
    void startSearch(int& N, int*& Rows) const;
    void finishSearch(int *const Rows, const QModelIndexList& Result) const;

    MTable *table;
    TableViewWindowCore *const mCore;
    int lastClickedHeaderIndex = -1;
};

#endif
