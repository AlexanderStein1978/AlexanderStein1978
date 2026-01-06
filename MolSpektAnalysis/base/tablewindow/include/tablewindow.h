//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLEWINDOW_H
#define TABLEWINDOW_H

#include "mdichild.h"
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
class ElState;


class TableWindow : public MDIChild
{
	Q_OBJECT
	
public:
	TableWindow(Type typ = TermEnergyTable, MainWindow *MW = 0, Molecule *M = 0);
    virtual ~TableWindow();
	virtual void setTabDimensions(int NRows, int NCols) = 0;
	virtual void getTabDimensions(int &NRows, int &NCols) = 0;
	virtual void setRowData(int Row, QString *Data) = 0;
	virtual void setMolecule(Molecule *Mol) = 0;
	Molecule *getMolecule();
	void setSource(QString nSource);
	QString getSource() const;
	int getvMax();
	int getJMax();
	void setvMax(int vM);
	void setJMax(int JM);
	double getError();
	bool isAssigned();
	virtual QString **getData(int &NRows, int &NCols) = 0;
	virtual QStringList getHorizontalHeaderLabels() = 0;
	virtual void setData(QString **Data, int NRows, int NCols) = 0;
	virtual void setHorizontalHeader(QStringList &Labels) = 0;
	virtual void setVerticalHeader(QStringList &Labels) = 0;
	virtual void cutRows(int &numRows, int &numColumns, QString **&Data) = 0;
	virtual void copyRows(int &numRows, int &numColums, QString **&Data) = 0;
	virtual void insertRows(int numRows, int numColumns, QString **Data) = 0;
	virtual void MarkLines(int *rN, int N) = 0;
	virtual void shiftCellValue(int n) = 0;
	virtual void exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns) = 0;
	void setViewnRows(MDIChild *Viewer, int NRows, int *Rows);
	virtual void search(int column, int value, int smeqla) = 0;
	virtual void search(int column, double value, int smeqla) = 0;
	virtual void search(QString Text, int column=-1, bool completeCell = false) = 0;
	virtual void setEditable(bool Editable) = 0;
	virtual void getViewnE(int *&Js, double *&E, int &N) = 0;
	virtual void DeleteRows() = 0;
	virtual void AddRow() = 0;
	virtual void RemoveDoubled() = 0;
    virtual void setCellText(QString Text) = 0;
	inline virtual int getNumLines() const = 0;
	inline virtual int getNumColumns() const = 0;
	
	inline virtual ElState *getElState()
	{
		return nullptr;
	}

public slots:
	void setName(QString nName);
	void setName();
	bool writeData(QString Filename = "");
	bool readData(QString Filename = "");
	virtual void tabItemChanged(QTableWidgetItem *Item) = 0;
	
protected:
	void resizeEvent(QResizeEvent *e);
	virtual void resizeHelper() = 0;
	virtual void setIsoIcon(int Col, int c = 0) = 0;
    virtual int *heapSort(bool sortFuncs(const QTableWidget *const Tab, const int n, const int m)) const = 0;
	virtual void sortTab(int *SortArray) = 0;
	void getViewnRows(bool *RV);
	void setNumIterations(int Finished, int Max = -1);
	void setNumParIt(int N);
	int getNumParIt();
	void setMaxParFits(int Max);
    virtual bool checkAllConnections(int FileColumn) = 0;
    virtual void shrinkAllSpectRefs(int FileColumn) = 0;
	virtual bool readData(QTextStream& S) = 0;
	virtual inline void writeData(QTextStream&) {}

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
	QList<ViewList> ViewLists;
	QProgressBar *Progress;
	Type Typ;

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
};

#endif
