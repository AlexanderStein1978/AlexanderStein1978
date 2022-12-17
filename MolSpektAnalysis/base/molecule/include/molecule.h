//
// C++ Interface: molecule
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MOLECULE_H
#define MOLECULE_H

#include <QComboBox>
#include <QLabel>
#include <QTableWidget>
#include <QPushButton>

#include "constants.h"

#include "elstate.h"
#include "MainWindow.h"
#include "mdichild.h"


class LineTable;
class IsoTab;

class Molecule : public MDIChild
{
	Q_OBJECT
public:
    Molecule(MainWindow *mw);
    ~Molecule();
	bool readData(QString Filename);
	Atom *getAtom1();
	Atom *getAtom2();
	void setAtoms(Atom *A1, Atom *A2);
	QString getState(int Index);
	ElState *getStateP(int Index);
	ElState *getState(QString Name);
    int addState(QString name, const bool i_reading);
	int addTransition(ElState *LState = 0, ElState *UState = 0, LineTable *Tab = 0, 
					  FCFTab *FCFT = 0);
	void addTerm(int SIndex, TermTable *TermTab);
	void addDunK(int SIndex, DunTable *DunTab);
	TermTable *getTerm(int SIndex);
    TermTable *getTermTable(QString Name);
    ElState *getStateOfTermTable(QString Name);
	DunTable *getDunK(int SIndex);
	void addPot(int SIndex, Potential *Pot);
	Potential *getPot(int SIndex);
	Potential *getPot(QString SName);
	QString getPotName(int SIndex);
	void addFitData(int SIndex, FitData *FD);
	FitData *getFitData(int SIndex);
    FitData *getFitData(QString Name);
	void addLines(int SIndex1, int SIndex2, LineTable *Lines);
	void addLineTable(LineTable *Tab, ElState *LState = 0, ElState *UState = 0);
	LineTable *getLines(int SIndex1, int SIndex2);
	LineTable *getLineTable(ElState *lState, ElState *uState);
	LineTable *getLineTable(QString Name);
	void addFCF(int SIndex1, int SIndex2, FCFTab *FCF);
	void addFCFTable(FCFTab *Tab, ElState *LState = 0, ElState *UState = 0);
	FCFTab *getFCF(int SIndex1, int SIndex2);
	FCFTab *getFCFTable(ElState *lState, ElState *uState);
	bool FCFavailable();
	int getNumStates();
	int getNumTransitions();
	double getTransitionStrength(ElState *State1, ElState *State2);
	Transition *getTransitionP(int Index);
	IsoTab *getIso();
	int getNumIso();
	int getJStep(int Iso);
	double getIsoMass(int Iso);
	void getRefIso(int &At1, int &At2);
	void setRefIso(int A1, int A2);
	void UpdateTermBox(QComboBox *B = 0);
	void UpdateDunBox(QComboBox *B = 0);
	void UpdatePotBox(QComboBox *B = 0);
	void UpdateFitDataBox(QComboBox *B = 0);
	void UpdateLineBox(QComboBox *B = 0);
	void UpdateFCFBox(QComboBox *B = 0);
	//void UpdateStateBox(QComboBox *B, QString *T = 0, QString C = "unknown");
	void UpdateStateBox();
	void ChangeObjName(MDIChild::Type type);
	void getKnownLevels(IsoTab *&Iso, int &NStates, int &mv, int &mJ, bool ****&Levels);
	void getTermData(int StateNum, int &numComp, int &numIso, int &numv, int &numJ, double ****&Data);
    void updateTransitionName(int transIndex, int lowerStateIndex, int upperStateIndex);
	void setStateName(int n, QString Name);
    bool checkAllConnections();
    void shrinkAllSpectRefs();
	
public slots:
	void updateAtoms();
	bool writeData(QString Filename = "");
	void TStatesClicked(int row, int column);
	void TTransitionsClicked(int row, int column);
private slots:
	void updateNI();
	void TStatesChanged(QTableWidgetItem *I);
private:
	int getTransitionIndex(int SIndex1, int SIndex2);
	Atom *atom1, *atom2;
	ElState *States[MaxStates];
	Transition *Transitions[MaxTransitions];
	QComboBox *BAtom1, *BAtom2, *RefIso;
	QLabel *Name;
	QString RIN;
	int numStates, numTransitions;
	QTableWidget *TStates, *TTransitions;
	MainWindow *MW;
	bool UNI;
};

#endif
