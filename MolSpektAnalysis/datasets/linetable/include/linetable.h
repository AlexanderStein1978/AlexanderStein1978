//
// C++ Interface: linetable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINETABLE_H
#define LINETABLE_H

#include "elstate.h"
#include "tablewindow.h"
#include "linetablesortfunctions.h"

#include <qmessagebox.h>

class Molecule;
class MainWindow;
class Marker;
class Progression;
class IsoTab;
class Spektrum;
class Progression;
class TableLine;
class TermEnergy;

struct vsOListElement;

class QComboBox;
class QPushButton;
class QStringList;
class QTableWidgetItem;
class QTableWidgetSelectionRange;
//class QList<QTableWidgetSelectionRange>;


enum TableCols {CPN, Cvs, CJs, Cvss, CJss, CF, CWN, Cerr, CIso, CFile, CSNR, CDev, CC,
				CFCF, CEUp, CEav, CEUma, CEdJ, CCalc, COmC};
#define TableNormCols 13


class LineTable : public TableWindow
{
	Q_OBJECT
	
public:
    LineTable(MainWindow *MW = 0, Molecule *M = 0, Transition *T = 0);
    ~LineTable();
	void setTransition(Transition *T);
	Transition *getTransition();
	int getAnzahlLinien();
    void getLines(int **Zuordnung, double *Energien, double *Uncertainties);
	void getLines(const QString &Filename, double **Lines, int *numLines);
	void getLines(TableLine *&L, int &N);
	void getSortedLines(TableLine *&L, int &N, int SortOrder = 0);
    void getgoodLines(int &N, TableLine *&L, int *mv = 0, int mJ = 0, bool SortFunction(const QTableWidget *const, const int, const int) = sortIJvP);
	void setUncertainty(int *RowNumbers, double *NewUncertainty, int NLines);
	int getNgL(int *mv = 0, int mJ = 0);
	int getNgTE(int *mv = 0, int mJ = 0);
	void getgoodTE(int &N, TermEnergy *&TE, int *mv = 0, int mJ = 0);
	int getMaxvs();
	int getMaxvss();
	int getMaxJs();
	int getMaxJss();
	Progression getSelectedProgression();
	void findSimilarProgression(Progression P);
	void getObsIso(bool *ObsIso, int NIso);
    void getProgressions(int &N, Progression *&P, bool SortFunction(const QTableWidget *const, const int, const int) = sortIJvP);
	bool writeData(QString Filename = "");
	bool readData(QString Filename = "");
	void addData(QString **Data, int NR, int NC);
	void setData(int NRows, int nCols, QStringList &HHeader, QTableWidgetItem ***Data);
	void splitTable();
	void WriteTFGS(QString FileName = "", int vso = 1000, vsOListElement *vsOList = 0);
	void W2AI(QString sldc = "");
	bool writeExcPotFitInput(QString FileName);
	void Assignvs();
	void AssignFC();
	void AddMarked(int AnzahlMarker, Marker *marker, Marker *LaserLine, QString SpektFile);
	void AcceptAssignments(QString SpektFile, bool accept);
    void MarkSelected();
	void MarkLines(int *Iso, int *vs, int *Js, double *WN, int N);
	void findErrors();
	
	inline void MarkLines(int *R, int N)
	{
		TableWindow::MarkLines(R, N);
	}
	
	inline ElState *getElState()
	{
		return (transition != 0 ? transition->getUpperState() : 0); 
	}
	
	void TakeOnChanges();
    void TestProgressions(int NumWFPoints);
    void DeleteRows();
	void deleteRows(int *rows, int N);
	void cutRows(int &numRows, int &numColumns, QString **&Data);
    bool ShowUpTerm();
    void ShowCalcRelInt(int NumWFPoints);
	void ShowGSDeviations();
    void FindBigDiff();
	void ShowWeakProgressions();
    TableWindow *ShowUpTermTable();	
	void Updatevs(int *nvs);
	void RemoveDoubled();
	void sortUpTermIvJ();
	void SortProg();
	void SortIJvP();
	void SortIvPJ();
    void SortFPInt();
	void SortSpectrum();
	void SortfRemDoubled();
	void SetError();
	void SetError(int NL, int *PN, int *vss, int *Jss, double *Err);
	void SetError(int NL, int *PN, int *vss, int *Jss, QString *Err);
	void DeleteRows(int NL, int *PN, int *vss = 0, int *Jss = 0);
	void SetPN();
	void Shiftvup();
	void Shiftvdown();
	void ShiftJup();
	void ShiftJdown();
	void Shiftvsup();
	void Shiftvsdown();
	void ShiftIso();
	void SetvssAscending();
	void Delete();
	void setvs();
	void setFC(QString nF);
	void sortbyvs();
	void getKnownLevels(int NI, int &mvs, int &mvss, int &mJs, int &mJss, bool ***&uL, bool ***&lL);
	void setMolecule(Molecule *mol);
	void getSelData(int *&Js, double *&E, int &N);
	void getViewnE(int *&Js, double *&E, int &N);

    inline bool checkAllConnections()
    {
        return TableWindow::checkAllConnections(CFile);
    }

    inline void shrinkAllSpectRefs()
    {
        TableWindow::shrinkAllSpectRefs(CFile);
    }

    inline void sortByProgNumber()
    {
        sortTab(heapSort(sortByProgression));
    }
	
public slots:
	void updateMarker(Spektrum *Spectrum);
    //void printf(const char* arg1);
	
private slots:
	void TabSelChanged();
	
signals:
	void DataChanged();
	
private:
	void sortTab(int *SortArray);
	void UpdateMarker(Spektrum *Spectrum = 0, int numLines = -1, int *Lines = 0, bool remove = false);
	
	Transition *transition;
	int NR, lRow, MaxPN, mvs, mJs, mIso, NpProg, NpL, *SelJs, NSel, *SO, NSO;
	double Error, OvError, *SelE;
	QString InFile, mSpectrum;
	QStringList HeaderLabels;
	TableWindow *termTable;
	IsoTab *Iso;
	QList<QTableWidgetSelectionRange> SelR;
};


#endif
