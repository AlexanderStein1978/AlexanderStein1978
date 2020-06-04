//
// C++ Interface: tablewindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLEWINDOW_H
#define TABLEWINDOW_H

#include "mdichild.h"
#include "MainWindow.h"
#include "viewlist.h"
#include "mtable.h"
#include "tabsortfunctor.h"

#include <QTableWidget>
#include <QString>
#include <QLineEdit>
#include <QLabel>
#include <QComboBox>
#include <QPushButton>
#include <QTableView>

class QProgressBar;
class Molecule;


class TableWindow : public MDIChild
{
	Q_OBJECT
	
public:
	TableWindow(Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    virtual ~TableWindow();
	void setTabDimensions(int NRows, int NCols);
	void getTabDimensions(int &NRows, int &NCols);
	void setRowData(int Row, QString *Data);
	virtual void setMolecule(Molecule *Mol);
	Molecule *getMolecule();	
	void setSource(QString nSource);
	QString getSource();
	int getvMax();
	int getJMax();
	void setvMax(int vM);
	void setJMax(int JM);
	double getError();
	bool isAssigned();
	QString **getData(int &NRows, int &NCols);
	QStringList getHorizontalHeaderLabels();
	void setData(QString **Data, int NRows, int NCols);
	void setHorizontalHeader(QStringList &Labels);
	void setVerticalHeader(QStringList &Labels);
	virtual void cutRows(int &numRows, int &numColumns, QString **&Data);
	void copyRows(int &numRows, int &numColums, QString **&Data);
	virtual void insertRows(int numRows, int numColumns, QString **Data);
	void MarkLines(int *rN, int N);
	void shiftCellValue(int n);
	void exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns);
	void setViewnRows(MDIChild *Viewer, int NRows, int *Rows);
	void search(int column, int value, int smeqla);
	void search(int column, double value, int smeqla);
	void search(QString Text, int column=-1, bool completeCell = false);
	void setEditable(bool Editable);
	virtual void getViewnE(int *&Js, double *&E, int &N);
	virtual void DeleteRows();
	virtual void AddRow();
	virtual void RemoveDoubled();
    void setCellText(QString Text);
	
	inline virtual ElState *getElState()
	{
		return 0;
	}
	
    inline virtual int getNumLines() const
	{
		if (Tab != 0) return Tab->rowCount();
		return 0;
	}
	
	inline int getNumColumns()
	{
		if (Tab != 0) return Tab->columnCount();
		return 0;
	}
	
public slots:
	void setName(QString nName);
	void setName();
	bool writeData(QString Filename = "");
	bool readData(QString Filename = "");
	void tabItemChanged(QTableWidgetItem *Item);
	
protected:
	void resizeEvent(QResizeEvent *e);
	void setIsoIcon(int Col, int c = 0);
    int *heapSort(bool sortFuncs(const QTableWidget *const Tab, const int n, const int m)) const;
	virtual void sortTab(int *SortArray);
	void getViewnRows(bool *RV);
	void setNumIterations(int Finished, int Max = -1);
	void setNumParIt(int N);
	int getNumParIt();
	void setMaxParFits(int Max);
    bool checkAllConnections(int FileColumn);
    void shrinkAllSpectRefs(int FileColumn);

    virtual inline QRegExp GetStartSpecialPartRegExp() const
    {
        return QRegExp();
    }

    virtual inline bool ReadSpecialPart(QTextStream& /*i_stream*/, const QString& /*i_startString*/)
    {
        return true;
    }
	
	Molecule *molecule;	
	QTableWidget *Tab;
	MTable *table;
	QLabel *vMLabel, *JMLabel, *errLabel, *IsoLabel, *CompLabel, *ViewLabel;
	QLabel *MolLabel, *USLabel, *LSLabel, *JsLabel, *JssLabel, *TLabel;
	QComboBox *Iso, *Comp, *View, *MolBox, *USBox, *LSBox, *TUnit, *JsB, *JssB;
	QLineEdit *vMax, *JMax, *error, *Js, *Jss, *Temp, *Date, *Pot1, *Pot2, *NumParFits;
	QPushButton *Calc;
	QString Filename;
	QPixmap getIsoIcon(int N);
	QList<ViewList> ViewLists;
	QProgressBar *Progress;

signals:
	void SourceChanged();
	void TabRowChanged(QString *Items, int N);
	void SelChanged();
	
private slots:
	void sourceChanged();
	void vMaxChanged();
	void JMaxChanged();
	void errorChanged();
	
private:
	QLineEdit *Source, *Name;
	QString tSource, tvMax, tJMax, terror;
	Type Typ;
	QPixmap *IsoIcon;
};

#endif
