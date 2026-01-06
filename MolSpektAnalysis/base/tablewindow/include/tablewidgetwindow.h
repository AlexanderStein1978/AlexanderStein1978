//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_H
#define TABLE_VIEW_WINDOW_H


#include "tablewindow.h"


class TableWidgetWindow : public TableWindow
{
public:
    TableWidgetWindow(Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    virtual ~TableWidgetWindow();

    void MarkLines(int *rN, int N) override;
    void setIsoIcon(int Col, int c = 0) override;
    void tabItemChanged(QTableWidgetItem *Item) override;
    void exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns) override;
    QString **getData(int &NRows, int &NCols) override;
    void setCellText(QString Text) override;
    void setData(QString **Data, int NRows, int NCols) override;
    void setEditable(bool Editable) override;
    void setHorizontalHeader(QStringList &Labels) override;
    void setVerticalHeader(QStringList &Labels) override;
    QStringList getHorizontalHeaderLabels() override;
    void search(int column, int value, int smeqla) override;
	void search(int column, double value, int smeqla) override;
	void search(QString Text, int column=-1, bool completeCell = false) override;

protected:
    void shiftCellValue(int n) override;
    int *heapSort(bool sortFuncs(const QTableWidget *const Tab, const int n, const int m)) const override;
    void shrinkAllSpectRefs(int FileColumn) override;
    void sortTab(int *SortArray) override;
};

#endif
