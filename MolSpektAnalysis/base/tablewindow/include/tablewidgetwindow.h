//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_WIDGET_WINDOW_H
#define TABLE_WIDGET_WINDOW_H


#include "tablewindow.h"


class TableWidgetWindow : public TableWindow
{
public:
    TableWidgetWindow(Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    virtual ~TableWidgetWindow();

    void setTabDimensions(int NRows, int NCols) override;
    void getTabDimensions(int &NRows, int &NCols) override;
    void setRowData(int Row, QString *Data) override;
    void AddRow() override;
    void DeleteRows() override;
    void cutRows(int &numRows, int &numColumns, QString **&Data) override;
	void copyRows(int &numRows, int &numColums, QString **&Data) override;
	void insertRows(int numRows, int numColumns, QString **Data) override;
    void MarkLines(int *rN, int N) override;
    void setIsoIcon(int Col, int c = 0);
    void tabItemChanged(QTableWidgetItem *Item);
    void exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns) override;
    QString **getData(int &NRows, int &NCols) override;
    void setCellText(QString Text);
    void setData(QString **Data, int NRows, int NCols) override;
    void setEditable(bool Editable) override;
    void setHorizontalHeader(QStringList &Labels) override;
    void setVerticalHeader(QStringList &Labels) override;
    QStringList getHorizontalHeaderLabels() override;
    void search(int column, int value, int smeqla) override;
	void search(int column, double value, int smeqla) override;
	void search(QString Text, int column=-1, bool completeCell = false) override;

    inline int getNumLines() const override
	{
		if (Tab != 0) return Tab->rowCount();
		return 0;
	}

	inline int getNumColumns() const override
	{
		if (Tab != 0) return Tab->columnCount();
		return 0;
	}

protected:
    void resizeHelper(QRect& G) override;
    bool readData(QTextStream& S) override;
    void writeData(QTextStream& S) override;
    void shiftCellValue(int n) override;
    int *heapSort(bool sortFuncs(const QTableWidget *const Tab, const int n, const int m)) const;
    void shrinkAllSpectRefs(int FileColumn);
    void sortTab(int *SortArray) override;
    bool checkAllConnections(int FileColumn);

    QTableWidget *Tab;
};

#endif
