//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINETABLE_H
#define LINETABLE_H


#include "Transition.h"
#include "tableviewwindow.h"
#include "linetablesortfunctions.h"
#include "linetablecore.h"

#include <qmessagebox.h>

class Molecule;
class MainWindow;
class Marker;
class Progression;
struct IsoTab;
class Spektrum;
class TableLine;
class TermEnergy;

struct vsOListElement;

class QComboBox;
class QPushButton;
class QStringList;
class QTableWidgetItem;
class QTableWidgetSelectionRange;
//class QList<QTableWidgetSelectionRange>;


class LineTable : public TableViewWindow
{
	Q_OBJECT
	
public:
    LineTable(MainWindow *MW = 0, Molecule *M = 0, Transition *T = 0);
    ~LineTable();
	void setTransition(Transition *T);
	Transition *getTransition();
	int getAnzahlLinien();
    void getLines(int **Zuordnung, double *Energien, double *Uncertainties);
	// void getLines(const QString &Filename, double **Lines, int *numLines);
	void getLines(TableLine *&L, int &N);
	void getSortedLines(TableLine *&L, int &N, int SortOrder = 0);
    void getgoodLines(int &N, TableLine *&L, int *mv = 0, int mJ = 0, bool SortFunction(const TableViewWindowCore *const, const int, const int) = sortIJvP);
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
    void getProgressions(int &N, Progression *&P, bool SortFunction(const TableViewWindowCore *const, const int, const int) = sortIJvP);
	bool writeData(QString Filename = "") override;
	bool readData(QString Filename = "") override;
	void addData(QString **Data, int NR, int NC);
	void setData(const std::vector<LineTableCore::TableCols>& mColumns, const std::vector<BaseData*> data);
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
	void setFC(const int nF);
	void sortbyvs();
	void getKnownLevels(int NI, int &mvs, int &mvss, int &mJs, int &mJss, bool ***&uL, bool ***&lL);
	void setMolecule(Molecule *mol);
	void getSelData(int *&Js, double *&E, int &N);
	void getViewnE(int *&Js, double *&E, int &N);
	void shrinkAllSpectRefs();

	inline void setData(const int row, BaseData* const data)
	{
		mCore->setData(row, data);
	}

	inline void MarkLines(int *R, int N) override
	{
		TableWindow::MarkLines(R, N);
	}

	inline ElState *getElState() override
	{
		return (transition != 0 ? transition->getUpperState() : 0);
	}

    inline bool checkAllConnections()
    {
        return TableViewWindow::checkAllConnections();
    }

    inline void sortByProgNumber()
    {
        sortTab(mCore->heapSort(sortByProgression));
    }

    inline QString getInfile() const
	{
		return InFile;
	}

	inline int getNumNewLines() const
	{
		return NpL;
	}
	
public slots:
	void updateMarker(Spektrum *Spectrum);
    //void printf(const char* arg1);
	
private slots:
	void TabSelChanged();
	void HeaderItemDoubleClicked(const int index);
	
signals:
	void DataChanged();
	
private:
	void sortTab(int *SortArray) override;
	void UpdateMarker(Spektrum *Spectrum = 0, int numLines = -1, int *Lines = 0, bool remove = false);
	
	Transition *transition;
	int NR, lRow, MaxPN, mvs, mJs, mIso, NpProg, NpL, *SelJs, NSel, *SO, NSO;
	double Error, OvError, *SelE;
	QString InFile, mSpectrum;
	TableWindow *termTable;
	IsoTab *Iso;
};


#endif
