//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QMouseEvent>
#include <QKeyEvent>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QPainter>
#include <QFont>
#include <QWidget>
#include <QAction>
#include <QMenu>
#include <QGridLayout>
#include <QPixmap>
#include <QPushButton>
#include <QTabWidget>
#include <QMessageBox>

#include <stdio.h>

#include "elstate.h"
#include "MainWindow.h"
#include "termtable.h"
#include "duntable.h"
#include "molecule.h"
#include "potential.h"
#include "fitdata.h"
#include "Combobox.h"
#include "Table.h"


ElState::ElState()
{
	stateNum = 0;
	termTable = 0;
	DunK = 0;
	Pot = 0;
	FitD = 0;
	MW = 0;
	TermB = 0;
	DunB = 0;
	PotB = 0;
	FDB = 0;
	TermTB = DunTB = PotTB = FDTB = 0;
	lambda = 0;
	mDun = 0;
	Mol = 0;
	mPot = 0;
	mTerm = 0;
	mFD = 0;
	numDun = 0;
	numPot = 0;
	numTerm = 0;
	numFD = 0;
	Omega = 0;
	Par = 0;
	S = 0;
    m_Be = 0.0;
	Sym = 0;
	HS = false;
	ShowTermMenu = ShowDunMenu = ShowPotMenu = ShowFitDataMenu = 0;
	ShowTermAction = ShowTMAction = ShowPotAction = ShowPMAction = ShowDunAction = ShowDMAction = 0;
	ShowFitDataAction = ShowFDMAction = 0;
	TBlock = false;
	QFont F;
	QBrush BC = palette().window();
	setWindowTitle("State properties:");
	setMinimumSize(300, 450);
	QGridLayout *Layout = new QGridLayout(this);
	QLabel *LL = new QLabel("Name:", this);
	Layout->addWidget(LL, 0, 0, 1, 1);
	Layout->setColumnMinimumWidth(0, 40);
	Layout->setColumnStretch(0, 2);
	StateName = new QLineEdit("A", this);
	Layout->addWidget(StateName, 0, 1, 1, 1);
	Layout->setColumnMinimumWidth(1, 40);
	Layout->setColumnStretch(1, 2);
	Layout->setColumnMinimumWidth(2, 20);
	Layout->setColumnStretch(2, 1);
	connect(StateName, SIGNAL(editingFinished()), this, SLOT(TypeChanged()));
	QString lambda = "\\lambda:";
	int w, h;
	QPixmap lambdaP(w = TextWidth(F, lambda), h = TextHeight(F, lambda));
	QPainter P(&lambdaP);
	P.fillRect(0, 0, w, h, BC);
	WriteText(P, 0, h, lambda, F, 0);
	QLabel *lambdaL = new QLabel(this);
	lambdaL->setPixmap(lambdaP);
	Layout->addWidget(lambdaL, 0, 3, 1, 1);
	Layout->setColumnMinimumWidth(3, 40);
	Layout->setColumnStretch(3, 2);
	lambdaB = new QComboBox(this);
	Layout->addWidget(lambdaB, 0, 4, 1, 1);
	Layout->setColumnMinimumWidth(4, 40);
	Layout->setColumnStretch(4, 2);
	lambdaB->addItem("0");
	lambdaB->addItem("1");
	lambdaB->addItem("2");
	lambdaB->setEditable(false);
	connect(lambdaB, SIGNAL(currentIndexChanged(int)), this, SLOT(TypeChanged()));
	Layout->setColumnMinimumWidth(5, 20);
	Layout->setColumnStretch(5, 1);
	QLabel *SL = new QLabel("S:", this);
	Layout->addWidget(SL, 0, 6, 1, 1);
	Layout->setColumnMinimumWidth(6, 40);
	Layout->setColumnStretch(6, 2);
	SB = new QComboBox(this);
	Layout->addWidget(SB, 0, 7, 1, 1);
	Layout->setColumnMinimumWidth(7, 40);
	Layout->setColumnStretch(7, 2);
	SB->setEditable(false);
	connect(SB, SIGNAL(currentIndexChanged(int)), this, SLOT(TypeChanged()));
	QLabel *PL = new QLabel("Symmetry:", this);
	Layout->addWidget(PL, 1, 0, 1, 1);
	ParB = new QComboBox(this);
	Layout->addWidget(ParB, 1, 1, 1, 1);
	ParB->addItem("");
	ParB->addItem("g");
	ParB->addItem("u");
	ParB->setEditable(false);
	connect(ParB, SIGNAL(currentIndexChanged(int)), this, SLOT(TypeChanged()));
	QLabel *SyL = new QLabel("Symmetry:", this);
	Layout->addWidget(SyL, 1, 3, 1, 1);
	SymB = new QComboBox(this);
	Layout->addWidget(SymB, 1, 4, 1, 1);
	SymB->addItem("");
	SymB->addItem("+");
	SymB->addItem("-");
	SymB->setEditable(false);
	connect(SymB, SIGNAL(currentIndexChanged(int)), this, SLOT(TypeChanged()));
	QString Omega = "\\Omega:";
	QPixmap OmegaP(w = TextWidth(F, Omega), h = TextHeight(F, Omega));
	QPainter OP(&OmegaP);
	OP.fillRect(0, 0, w, h, BC);
	WriteText(OP, 0, h, Omega, F, 0);
	QLabel *OmegaL = new QLabel(this);
	OmegaL->setPixmap(OmegaP);
	Layout->addWidget(OmegaL, 1, 6, 1, 1);
	OmegaB = new QComboBox(this);
	Layout->addWidget(OmegaB, 1, 7, 1, 1);
	OmegaB->addItem("0");
	OmegaB->addItem("1");
	OmegaB->addItem("2");
	OmegaB->setEditable(false);
	connect(OmegaB, SIGNAL(currentIndexChanged(int)), this, SLOT(TypeChanged()));
	Layout->setRowMinimumHeight(2, 20);
	Register = new QTabWidget(this);
	Layout->addWidget(Register, 3, 0, 1, 8);
	tTable = new Table(nTermR = 10, 3, this);
	QStringList HeaderLabels;
	tTable->setHorizontalHeaderLabels(HeaderLabels << "Name" << "File name" << "Source");
	connect(tTable, SIGNAL(cellChanged(int, int)), this, SLOT(termTableChanged(int, int)));
	connect(tTable, SIGNAL(cellClicked(int, int)), this, SLOT(termTableClicked(int, int)));
	connect(tTable, SIGNAL(RowDeleted(bool*)), this, SLOT(deleteTermTables(bool*)));
	Register->addTab(tTable, "Term energy sets");
	dTable = new Table(nDunR = 10, 3, this);
	dTable->setHorizontalHeaderLabels(HeaderLabels);
	connect(dTable, SIGNAL(cellChanged(int, int)), this, SLOT(dunTableChanged(int, int)));
	connect(dTable, SIGNAL(cellClicked(int, int)), this, SLOT(dunTableClicked(int, int)));
	connect(dTable, SIGNAL(RowDeleted(bool*)), this, SLOT(deleteDunTables(bool*)));
	Register->addTab(dTable, "Dunham coefficient sets");
	pTable = new Table(nPotR = 10, 3, this);
	pTable->setHorizontalHeaderLabels(HeaderLabels);
	connect(pTable, SIGNAL(cellChanged(int, int)), this, SLOT(potTableChanged(int, int)));
	connect(pTable, SIGNAL(cellClicked(int, int)), this, SLOT(potTableClicked(int, int)));
	connect(pTable, SIGNAL(RowDeleted(bool*)), this, SLOT(deletePotentials(bool*)));
	Register->addTab(pTable, "Potentials");
	fDTable = new Table(nFDR = 10, 3, this);
	fDTable->setHorizontalHeaderLabels(HeaderLabels);
	connect(fDTable, SIGNAL(cellChanged(int,int)), this, SLOT(fitDataTableChanged(int,int)));
	connect(fDTable, SIGNAL(cellClicked(int,int)), this, SLOT(fitDataTableClicked(int,int)));
	connect(fDTable, SIGNAL(RowDeleted(bool*)), this, SLOT(deleteFitDataSets(bool*)));
	Register->addTab(fDTable, "Fit datasets");
}

ElState::~ElState()
{
	//printf("ElState::~ElState\n");
	if (DunK != 0) delete[] DunK;
	if (termTable != 0) delete[] termTable;
	if (Pot != 0) delete[] Pot;
	if (FitD != 0) delete[] FitD;
	if (DunTB != 0) delete[] DunTB;
	if (PotTB != 0) delete[] PotTB;
	if (TermTB != 0) delete[] TermTB;
	if (FDTB != 0) delete[] FDTB;
	//printf("Ende ~ElState\n"); 
}

void ElState::addDunTable(DunTable *Table)
{
	addDunTable(Table, -1);
}

void ElState::addDunTable(DunTable *Table, int Index)
{
	QString cName, nName;
	int n;
	if (DunB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 2);
	TBlock = true;
	if (Index >= 0 && Index < numDun) cName = DunTB[Index]->currentText();
	if (DunK == 0) 
	{
		DunK = new DunTable*[10];
		for (n=0; n < 10; n++) DunK[n] = 0;
		emit Added(stateNum, 2);
	}
	if (DunTB == 0) DunTB = new ComboBox*[10];
	if (Index > numDun) Index = numDun;
	if (Index == nDunR)
	{
		DunTable **nDunK = new DunTable*[nDunR + 10];
		ComboBox **nDunTB = new ComboBox*[nDunR + 10];
		for (n=0; n < nDunR; n++)
		{
			nDunK[n] = DunK[n];
			nDunTB[n] = DunTB[n];
		}
		while (n < nDunR + 10) nDunK[n++] = 0;
		delete[] DunK;
		delete[] DunTB;
		DunK = nDunK;
		DunTB = nDunTB;
		dTable->setRowCount(nDunR += 10);
	}
	if (Index < 0) 
	{
		Index = numDun;
		dTable->setCellWidget(numDun, 0, DunTB[numDun] = new ComboBox(0));
		if (Table != 0)
		{
			DunTB[numDun]->addItem(Table->getName());
			DunTB[numDun]->setCurrentIndex(DunTB[numDun]->count() - 1);
		}
		connect(DunTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkDunTables()));
	}
	else if (DunK[Index] != 0)
	{
		disconnect(DunK[Index], 0, this, SLOT(checkDunTable()));
		disconnect(DunTB[Index], SIGNAL(DoubleClicked(int)), DunK[Index], SLOT(show()));
		DunK[Index]->setMolecule(0);
	}
	else disconnect(DunTB[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadDunTable(int)));
	DunK[Index] = Table;
	if (Table != 0)
	{
		DunTB[Index]->addItem(nName = Table->getName());
		DunTB[Index]->setCurrentIndex(DunTB[Index]->count() - 1);
		dTable->setItem(Index, 1, new QTableWidgetItem(Table->getFileName()));
		dTable->setItem(Index, 2, new QTableWidgetItem(Table->getSource()));
		Table->setElState(this);
		Table->setMolecule(Mol);
		connect(DunTB[Index], SIGNAL(DoubleClicked(int)), Table, SLOT(show()));
		connect(Table, SIGNAL(fileNameChanged()), this, SLOT(checkDunTables()));
		connect(Table, SIGNAL(nameChanged(QString)), this, SLOT(DunNameChanged()));
		connect(Table, SIGNAL(SourceChanged()), this, SLOT(checkDunTables()));
	}
	if (Index == numDun) numDun++;
	refreshDunBox();
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void ElState::addDunTable(QString Name, QString FileName, QString Source)
{
	if (DunB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 2);
	int n;
	TBlock = true;
	if (DunK == 0) 
	{
		DunK = new DunTable*[10];
		for (n=0; n < 10; n++) DunK[n] = 0;
		emit Added(stateNum, 2);
	}
	if (DunTB == 0) DunTB = new ComboBox*[10];
	if (numDun == nDunR)
	{
		DunTable **nDunK = new DunTable*[nDunR + 10];
		ComboBox **nDunTB = new ComboBox*[nDunR + 10];
		for (n=0; n < nDunR; n++)
		{
			nDunK[n] = DunK[n];
			nDunTB[n] = DunTB[n];
		}
		while (n < nDunR + 10) nDunK[n++] = 0;
		delete[] DunK;
		delete[] DunTB;
		DunK = nDunK;
		DunTB = nDunTB;
		dTable->setRowCount(nDunR += 10);
	}
	dTable->setCellWidget(numDun, 0, DunTB[numDun] = new ComboBox(0));
	DunTB[numDun]->addItem(Name);
	//printf("DunTB[%d]=%d\n", numDun, DunTB[numDun]);
	//printf("DunTB[%d]=%s\n", numDun, DunTB[numDun]->currentText().ascii());
	dTable->setItem(numDun, 1, new QTableWidgetItem(FileName));
	dTable->setItem(numDun, 2, new QTableWidgetItem(Source));
	connect(DunTB[numDun], SIGNAL(currentIndexChanged(int)), this, SLOT(checkDunTables()));
	connect(DunTB[numDun], SIGNAL(DoubleClicked(int)), this, SLOT(loadDunTable(int)));
	DunTB[numDun]->setNumber(numDun);
	DunK[numDun++] = 0;
	refreshDunBox();
	TBlock = false;
}

void ElState::addFitData(FitData* set)
{
	addFitData(set, -1);
}

void ElState::addFitData(FitData* set, int Index)
{
	QString cName, nName;
	int n;
	if (FDB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 4);
	TBlock = true;
	if (Index >= 0 && Index < numFD) cName = FDTB[Index]->currentText();
	if (FitD == 0) 
	{
		FitD = new FitData*[10];
		for (n=0; n < 10; n++) FitD[n] = 0;
		emit Added(stateNum, 4);
	}
	if (FDTB == 0) FDTB = new ComboBox*[10];
	if (Index > numFD) Index = numFD;
	if (numFD == nFDR)
	{
		FitData **nFitD = new FitData*[nFDR + 10];
		ComboBox **nFDTB = new ComboBox*[nFDR + 10];
		for (n=0; n < nFDR; n++)
		{
			nFitD[n] = FitD[n];
			nFDTB[n] = FDTB[n];
		}
		while (n < nFDR + 10) nFitD[n++] = 0;
		delete[] FitD;
		delete[] FDTB;
		FitD = nFitD;
		FDTB = nFDTB;
		fDTable->setRowCount(nFDR += 10);
	}
	if (Index < 0) 
	{
		Index = numFD;
		fDTable->setCellWidget(numFD, 0, FDTB[numFD] = new ComboBox(0));
		if (set != 0)
		{
			FDTB[numFD]->addItem(set->getName());
			FDTB[numFD]->setCurrentIndex(FDTB[numFD]->count() - 1);
		}
		connect(FDTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFitDataSets()));
	}
	else if (FitD[Index] != 0)
	{
		disconnect(FitD[Index], 0, this, SLOT(checkFitDataSets()));
		disconnect(FDTB[Index], SIGNAL(DoubleClicked(int)), FitD[Index], SLOT(show()));
		FitD[Index]->setMolecule(0);
	}
	else disconnect(FDTB[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadFitData(int)));
	FitD[Index] = set;
	if (set != 0)
	{
		FDTB[Index]->addItem(nName = set->getName());
		FDTB[Index]->setCurrentIndex(FDTB[Index]->count() - 1);
		fDTable->setItem(Index, 1, new QTableWidgetItem(set->getFileName()));
		fDTable->setItem(Index, 2, new QTableWidgetItem(set->getSource()));
        set->setMolecule(Mol);
        set->setElState(this);
		connect(FDTB[Index], SIGNAL(DoubleClicked(int)), set, SLOT(show()));
		connect(set, SIGNAL(fileNameChanged()), this, SLOT(checkFitDataSets()));
		connect(set, SIGNAL(nameChanged(QString)), this, SLOT(FitDataNameChanged()));
		connect(set, SIGNAL(SourceChanged()), this, SLOT(checkFitDataSets()));
	}
	if (Index == numFD) numFD++;
	refreshFitDataBox();
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void ElState::addFitData(QString Name, QString FileName, QString Source)
{
	if (FDB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 4);
	TBlock = true;
	int n;
	if (FitD == 0) 
	{
		FitD = new FitData*[10];
		for (n=0; n < 10; n++) FitD[n] = 0;
		emit Added(stateNum, 4);
	}
	if (FDTB == 0) FDTB = new ComboBox*[10];
	if (numFD == nFDR)
	{
		FitData **nFitD = new FitData*[nFDR + 10];
		ComboBox **nFDTB = new ComboBox*[nFDR + 10];
		for (n=0; n < nFDR; n++)
		{
			nFitD[n] = FitD[n];
			nFDTB[n] = FDTB[n];
		}
		while (n < nFDR + 10) nFitD[n++] = 0;
		delete[] FitD;
		delete[] FDTB;
		FitD = nFitD;
		FDTB = nFDTB;
		fDTable->setRowCount(nFDR += 10);
	}
	fDTable->setCellWidget(numFD, 0, FDTB[numFD] = new ComboBox(0));
	FDTB[numFD]->addItem(Name);
	fDTable->setItem(numFD, 1, new QTableWidgetItem(FileName));
	fDTable->setItem(numFD, 2, new QTableWidgetItem(Source));
	connect(FDTB[numFD], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFitDataSets()));
	connect(FDTB[numFD], SIGNAL(DoubleClicked(int)), this, SLOT(loadFitData(int)));
	FDTB[numFD]->setNumber(numFD);
	FitD[numFD++] = 0;
	refreshFitDataBox();
	TBlock = false;
}

void ElState::addPotential(Potential *pot)
{
	addPotential(pot, -1);
}

void ElState::addPotential(Potential *pot, int Index)
{
	//printf("ElState::addPotential: PotB=%d, Mol=%d\n", PotB, Mol);
	QString cName, nName;
	int n;
	if (PotB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 3); 
	TBlock = true;
	if (Index >= 0 && Index < numPot) cName = PotTB[Index]->currentText();
	if (Pot == 0) 
	{
		Pot = new Potential*[10];
		for (n=0; n < 10; n++) Pot[n] = 0;
		emit Added(stateNum, 3);
	}
	if (PotTB == 0) PotTB = new ComboBox*[10];
	if (Index > numPot) Index = numPot;
	if (numPot == nPotR)
	{
		Potential **nPot = new Potential*[nPotR + 10];
		ComboBox **nPotTB = new ComboBox*[nPotR + 10];
		for (n=0; n < nPotR; n++)
		{
			nPot[n] = Pot[n];
			nPotTB[n] = PotTB[n];
		}
		while (n < nPotR + 10) nPot[n++] = 0;
		delete[] Pot;
		delete[] PotTB;
		Pot = nPot;
		PotTB = nPotTB;
		pTable->setRowCount(nPotR += 10);
	}
	if (Index < 0)
	{
		Index = numPot;
		pTable->setCellWidget(numPot, 0, PotTB[numPot] = new ComboBox(0));
		if (pot != 0)
		{
			PotTB[numPot]->addItem(pot->getName());
			PotTB[numPot]->setCurrentIndex(PotTB[numPot]->count() - 1);
		}
		connect(PotTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkPotentials()));
	}
	else if (Pot[Index] != 0)
	{
		disconnect(Pot[Index], 0, this, SLOT(checkPotentials()));
		disconnect(PotTB[Index], SIGNAL(DoubleClicked(int)), Pot[Index], SLOT(show()));
		Pot[Index]->setMolecule(0);
	}
	else disconnect(Pot[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadPotential(int)));
	Pot[Index] = pot;
	if (pot != 0)
	{
		PotTB[Index]->addItem(nName = pot->getName());
		PotTB[Index]->setCurrentIndex(PotTB[Index]->count() - 1);
		pTable->setItem(Index, 1, new QTableWidgetItem(pot->getFileName()));
		pTable->setItem(Index, 2, new QTableWidgetItem(pot->getSource()));
		pot->setMolecule(Mol);
		pot->setElState(this);
		connect(PotTB[Index], SIGNAL(DoubleClicked(int)), pot, SLOT(show()));
		connect(pot, SIGNAL(fileNameChanged()), this, SLOT(checkPotentials()));
		connect(pot, SIGNAL(nameChanged(QString)), this, SLOT(PotNameChanged()));
		connect(pot, SIGNAL(SourceChanged()), this, SLOT(checkPotentials()));
	}
	if (Index == numPot) numPot++;
	refreshPotBox();
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void ElState::addPotential(QString Name, QString FileName, QString Source)
{
	if (PotB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 3);
	TBlock = true;
	int n;
	if (Pot == 0) 
	{
		Pot = new Potential*[10];
		for (n=0; n < 10; n++) Pot[n] = 0;
		emit Added(stateNum, 3);
	}
	if (PotTB == 0) PotTB = new ComboBox*[10];
	if (numPot == nPotR)
	{
		Potential **nPot = new Potential*[nPotR + 10];
		ComboBox **nPotTB = new ComboBox*[nPotR + 10];
		for (n=0; n < nPotR; n++)
		{
			nPot[n] = Pot[n];
			nPotTB[n] = PotTB[n];
		}
		while (n < nPotR + 10) nPot[n++] = 0;
		delete[] Pot;
		delete[] PotTB;
		Pot = nPot;
		PotTB = nPotTB;
		pTable->setRowCount(nPotR += 10);
	}
	pTable->setCellWidget(numPot, 0, PotTB[numPot] = new ComboBox(0));
	PotTB[numPot]->addItem(Name);
	connect(PotTB[numPot], SIGNAL(currentIndexChanged(int)), this, SLOT(checkPotentials()));
	connect(PotTB[numPot], SIGNAL(DoubleClicked(int)), this, SLOT(loadPotential(int)));
	PotTB[numPot]->setNumber(numPot);
	pTable->setItem(numPot, 1, new QTableWidgetItem(FileName));
	pTable->setItem(numPot, 2, new QTableWidgetItem(Source));
	Pot[numPot++] = 0;
	refreshPotBox();
	TBlock = false;
}

void ElState::addTermTable(TermTable *table)
{
	addTermTable(table, -1);
}

void ElState::addTermTable(TermTable *table, int Index)
{
	QString nName, cName;
	int n;
	if (TermB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 1);
	TBlock = true;
	//printf("ElState::addTermTable()\n");
	if (Index >= 0 && Index < numTerm) cName = TermTB[Index]->currentText();
	if (termTable == 0) 
	{
		termTable = new TermTable*[10];
		for (n=0; n < 10; n++) termTable[n] = 0;
		emit Added(stateNum, 1);
	}
	if (TermTB == 0) TermTB = new ComboBox*[10];
	if (Index > numTerm) Index = numTerm;
	if (numTerm == nTermR)
	{
		TermTable **ntermTable = new TermTable*[nTermR + 10];
		ComboBox **nTermTB = new ComboBox*[nTermR + 10];
		for (n=0; n < nTermR; n++)
		{
			ntermTable[n] = termTable[n];
			nTermTB[n] = TermTB[n];
		}
		while (n < nTermR + 10) ntermTable[n++] = 0;
		delete[] termTable;
		delete[] TermTB;
		termTable = ntermTable;
		TermTB = nTermTB;
		tTable->setRowCount(nTermR += 10);
	}
	if (Index < 0)
	{
		Index = numTerm++;
		tTable->setCellWidget(Index, 0, TermTB[Index] = new ComboBox(0));
		if (table != 0)
		{
			TermTB[Index]->addItem(table->getName());
			TermTB[Index]->setCurrentIndex(TermTB[Index]->count() - 1);
		}
		//printf("Nach setCurrentIndex\n");
		connect(TermTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkTermTables()));
	}
	else if (termTable[Index] != 0)
	{
		disconnect(termTable[Index], 0, this, SLOT(checkTermTable()));
		disconnect(TermTB[Index], SIGNAL(DoubleClicked(int)), termTable[Index], SLOT(show));
		termTable[Index]->setMolecule(0, 0);
	}
	else disconnect(TermTB[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadTermTable(int)));
	termTable[Index] = table;
	if (table != 0)
	{
		TermTB[Index]->addItem(nName = table->getName());
		TermTB[Index]->setCurrentIndex(TermTB[Index]->count() - 1);
		tTable->setItem(Index, 1, new QTableWidgetItem(table->getFileName()));
		tTable->setItem(Index, 2, new QTableWidgetItem(table->getSource()));
		table->setMolecule(Mol, this);
		//printf("nach table\n");
		connect(TermTB[Index], SIGNAL(DoubleClicked(int)), table, SLOT(show()));
		connect(table, SIGNAL(fileNameChanged()), this, SLOT(checkTermTables()));
		connect(table, SIGNAL(nameChanged(QString)), this, SLOT(TermNameChanged()));
		connect(table, SIGNAL(SourceChanged()), this, SLOT(checkTermTables()));
	}
	if (Index == numTerm) numTerm++;
	//printf("vor refreshTermBox\n");
	refreshTermBox();
	//printf("Ende addTermTable\n");
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void ElState::addTermTable(QString Name, QString FileName, QString Source)
{
	//printf("ElState::addTermTable(Name)\n");
	if (TermB == 0 && Mol != 0) Mol->TStatesClicked(stateNum, 1);
	TBlock = true;
	int n;
	if (termTable == 0) 
	{
		termTable = new TermTable*[10];
		for (n=0; n < 10; n++) termTable[n] = 0;
		emit Added(stateNum, 1);
	}
	if (TermTB == 0) TermTB = new ComboBox*[10];
	if (numTerm == nTermR)
	{
		TermTable **ntermTable = new TermTable*[nTermR + 10];
		ComboBox **nTermTB = new ComboBox*[nTermR + 10];
		for (n=0; n < nTermR; n++)
		{
			ntermTable[n] = termTable[n];
			nTermTB[n] = TermTB[n];
		}
		while (n < nTermR + 10) ntermTable[n++] = 0;
		delete[] termTable;
		delete[] TermTB;
		termTable = ntermTable;
		TermTB = nTermTB;
		tTable->setRowCount(nTermR += 10);
	}
	tTable->setCellWidget(numTerm, 0, TermTB[numTerm] = new ComboBox(0));
	TermTB[numTerm]->addItem(Name);
	connect(TermTB[numTerm], SIGNAL(currentIndexChanged(int)), this, SLOT(checkTermTables()));
	connect(TermTB[numTerm], SIGNAL(DoubleClicked(int)), this, SLOT(loadTermTable(int)));
	TermTB[numTerm]->setNumber(numTerm);
	tTable->setItem(numTerm, 1, new QTableWidgetItem(FileName));
	tTable->setItem(numTerm, 2, new QTableWidgetItem(Source));
	termTable[numTerm++] = 0;
	refreshTermBox();
	TBlock = false;
	//printf("Ende addTermTable\n");
}

void ElState::checkDunTables()
{
	if (MW == 0 || TBlock) return;
	int i, d, D = MW->getNumDunTables();
	bool Changed = false, TC = false;
	QString N, Text;
	DunTable *DT = 0;
	QTableWidgetItem *I;
	for (i=0; i < numDun; i++)
	{
		N = DunTB[i]->currentText();
		if (N.isEmpty())
		{
			bool Del[numDun];
			for (d=0; d < numDun; d++) Del[d] = true;
			Del[i] = false;
			deleteDunTables(Del);
			return;
		}
		else
		{
			if (DunK[i] != 0 ? N != DunK[i]->getName() : DunTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (DT = MW->getDunTable(d))->getName() != N : false); d++) ;
				if (d<D) 
				{
					addDunTable(DT, i);
					Changed = true;
					TC = true;
				}
			}
			if (DunK[i] != 0)
			{
				if ((I = dTable->item(i, 1))->text() != (Text = DunK[i]->getFileName())) 
				{
					I->setText(Text);
					TC = true;
				}
				if ((I = dTable->item(i, 2))->text() != (Text = DunK[i]->getSource())) 
				{
					I->setText(Text);
					TC = true;
				}
			}
		}
	}
	if (Changed) refreshDunBox();
	if (TC && Mol != 0) Mol->Changed();
}

void ElState::checkFitDataSets()
{
	if (MW == 0 || TBlock) return;
	int i, d, D = MW->getNumFitDataSets();
	bool Changed = false, TC = false;
	QString N, Text;
	FitData *FD = 0;
	QTableWidgetItem *I;
	for (i=0; i < numFD; i++)
	{
		N = FDTB[i]->currentText();
		if (N.isEmpty())
		{
			bool Del[numFD];
			for (d=0; d < numFD; d++) Del[d] = true;
			Del[i] = false;
			deleteFitDataSets(Del);
			return;
		}
		else
		{
			if (FitD[i] != 0 ? N != FitD[i]->getName() : FDTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (FD = MW->getFitData(d))->getName() != N : false); d++) ;
				if (d<D)
				{
					addFitData(FD, i);
					Changed = true;
					TC = true;
				}
			}
			if (FitD[i] != 0)
			{
				if ((I = fDTable->item(i, 1))->text() != (Text = FitD[i]->getFileName()))
				{
					I->setText(Text);
					TC = true;
				}
				if ((I = fDTable->item(i, 2))->text() != (Text = FitD[i]->getSource()))
				{
					I->setText(Text);
					TC = true;
				}
			}
		}
	}
	if (Changed) refreshFitDataBox();
	if (TC && Mol != 0) Mol->Changed();
}

void ElState::checkPotentials()
{
	if (MW == 0 || TBlock) return;
	//printf("ElState::checkPotentials\n");
	int i, d, D = MW->getNumPotentials();
	bool Changed = false, TC = false;
	QString N, Text;
	Potential *P = 0;
	QTableWidgetItem *I;
	for (i=0; i < numPot; i++)
	{
		N = PotTB[i]->currentText();
		if (N.isEmpty())
		{
			bool Del[numPot];
			for (d=0; d < numPot; d++) Del[d] = true;
			Del[i] = false;
			deletePotentials(Del);
			return;
		}
		else
		{
			if (Pot[i] != 0 ? N != Pot[i]->getName() : PotTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (P = MW->getPotential(d))->getName() != N : false); d++) ;
				if (d<D)
				{
					addPotential(P, i);
					TC = Changed = true;
				}
			}
			if (Pot[i] != 0)
			{
				//printf("vor FileName\n");
				if ((I = pTable->item(i, 1))->text() != (Text = Pot[i]->getFileName())) 
				{
					I->setText(Text);
					TC = true;
				}
				//printf("Nach Filename\n");
				if ((I = pTable->item(i, 2))->text() != (Text = Pot[i]->getSource())) 
				{
					I->setText(Text);
					TC = true;
				}
				//printf("Nach Source\n");
			}
		}
	}
	if (Changed) refreshPotBox();
	if (TC && Mol != 0) Mol->Changed();
	//printf("Ende checkPotentials\n");
}

void ElState::checkTermTables()
{
	//printf("ElState::checkTermTables()\n");
	if (MW == 0 || TBlock) return;
	int i, d, D = MW->getNumTermTables();
	bool Changed = false, TC = false;
	QString N, Text;
	TermTable *T = 0;
	QTableWidgetItem *I;
	//printf("D=%d\n", D);
	for (i=0; i < numTerm; i++)
	{
		N = TermTB[i]->currentText();
		//printf("i=%d, numTerm=%d, N=%s\n", i, numTerm, N.ascii());
		if (N.isEmpty())
		{
			bool Del[numTerm];
			for (d=0; d < numTerm; d++) Del[d] = true;
			Del[i] = false;
			deleteTermTables(Del);
			return;
		}
		else
		{
			if (termTable[i] != 0 ? N != termTable[i]->getName() : TermTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (T = MW->getTermTable(d))->getName() != N : false); d++) ;
				//printf("TermTable[0]=%s, N=%s, d=%d\n", 
					//   MW->getTermTable(0)->getName().ascii(), N.ascii(), d);
				if (d<D)
				{
					//printf("d=%d, D=%d, Text=%s, Name=%s\n", 
						//   d, D, N.toAscii().data(), T->getName().toAscii().data());
					addTermTable(T, i);
					Changed = TC = true;
				}
			}
			if (termTable[i] != 0)
			{
				if ((I = tTable->item(i, 1))->text() != (Text = termTable[i]->getFileName()))
				{
					I->setText(Text);
					TC = true;
				}
				if ((I = tTable->item(i, 2))->text() != (Text = termTable[i]->getSource()))
				{
					I->setText(Text);
					TC = true;
				}
			}
		}
	}
	if (Changed) refreshTermBox();
	if (TC && Mol != 0) Mol->Changed();
	//printf("Ende checkTermTables\n");
}

void ElState::deleteDunTables(bool *delR)
{
	int i, j, mD = mDun;
	TBlock = true;
	for (i=0; i < numDun; i++) if (!delR[i])
	{
		if (i == mDun) mDun = 0;
		if (DunK[i] == 0) 
			disconnect(DunTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadDunTable(int)));
		else
		{
			disconnect(DunK[i], 0, this, SLOT(checkDunTables()));
			disconnect(DunTB[i], SIGNAL(DoubleClicked(int)), DunK[i], SLOT(show()));
			DunK[i]->setMolecule(0);
		}
		DunK[i] = 0;
		if (dTable->item(i, 1) != 0) dTable->item(i, 1)->setText("");
		if (dTable->item(i, 2) != 0) dTable->item(i, 2)->setText("");
	}
	for (i=0; (i < numDun ? delR[i] : false); i++) ;
	for (j=i+1; j < numDun; j++) if (delR[j]) 
	{
		if (mDun == j) mDun = i;
		DunTB[i]->insertItem(DunTB[i]->count() - 1, DunTB[j]->currentText());
		DunTB[j]->setCurrentIndex(DunTB[j]->count() - 1);
		DunTB[i]->setCurrentIndex(DunTB[i]->count() - 2);
		if (DunK[j] != 0) 
		{
			disconnect(DunTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(DunTB[i], SIGNAL(DoubleClicked(int)), DunK[j], SLOT(show()));
		}
		else if (DunK[i] != 0)
		{
			disconnect(DunTB[i], SIGNAL(DoubleClicked(int)), DunK[i], SLOT(show()));
			connect(DunTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadDunTable(int)));
		}
		dTable->setItem(i, 1, dTable->takeItem(j, 1));
		dTable->setItem(i, 2, dTable->takeItem(j, 2));
		DunK[i++] = DunK[j];
	}
	TBlock = false;
	numDun = i;
	refreshDunBox();
	if (mD != mDun) 
	{
		emit DunhamChanged();
		if (DunB != 0) DunB->setCurrentIndex(mDun);
	}
	if (Mol != 0) Mol->Changed();
}

void ElState::deleteFitDataSets(bool* delR)
{
	int i, j, mF = mFD;
	TBlock = true;
	for (i=0; i < numFD; i++) if (!delR[i])
	{
		if (i == mFD) mFD = 0;
		if (FitD[i] == 0) 
			disconnect(FDTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadFitData(int)));
		else
		{
			disconnect(FitD[i], 0, this, SLOT(checkFitDataSets()));
			disconnect(FDTB[i], SIGNAL(DoubleClicked(int)), FitD[i], SLOT(show()));
			FitD[i]->setMolecule(0);
		}
		FitD[i] = 0;
		if (fDTable->item(i, 1) != 0) fDTable->item(i, 1)->setText("");
		if (fDTable->item(i, 2) != 0) fDTable->item(i, 2)->setText("");
	}
	for (i=0; (i < numFD ? delR[i] : false); i++) ;
	for (j=i+1; j < numFD; j++) if (delR[j])
	{
		if (mFD == j) mFD = i;
		FDTB[i]->insertItem(FDTB[i]->count() - 1, FDTB[j]->currentText());
		FDTB[i]->setCurrentIndex(FDTB[i]->count() - 2);
		FDTB[j]->setCurrentIndex(FDTB[j]->count() - 1);
		if (FitD[j] != 0)
		{
			disconnect(FDTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(FDTB[i], SIGNAL(DoubleClicked(int)), DunK[j], SLOT(show()));
		}
		else if (FitD[i] != 0)
		{
			disconnect(FDTB[i], SIGNAL(DoubleClicked(int)), FitD[i], SLOT(show()));
			connect(FDTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadFitData(int)));
		}
		fDTable->setItem(i, 1, fDTable->takeItem(j, 1));
		fDTable->setItem(i, 2, fDTable->takeItem(j, 2));
		FitD[i++] = FitD[j];
	}
	TBlock = false;
	numFD = i;
	if (mF != mFD)
	{
		emit FitDataChanged();
		if (FDB != 0) FDB->setCurrentIndex(mFD);
	}
	if (Mol != 0) Mol->Changed();
}

void ElState::deletePotentials(bool *delR)
{
	int i, j, mP = mPot;
	TBlock = true;
	for (i=0; i < numPot; i++) if (!delR[i])
	{
		if (mPot == i) mPot = 0;
		if (Pot[i] == 0) 
			disconnect(PotTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadPotential(int)));
		else
		{
			disconnect(Pot[i], 0, this, SLOT(checkPotentials()));
			disconnect(PotTB[i], SIGNAL(DoubleClicked(int)), Pot[i], SLOT(show()));
			Pot[i]->setMolecule(0);
		}
		Pot[i] = 0;
		if (pTable->item(i, 1) != 0) pTable->item(i, 1)->setText("");
		if (pTable->item(i, 2) != 0) pTable->item(i, 2)->setText("");
	}
	for (i=0; (i < numPot ? delR[i] : false); i++) ;
	for (j=i+1; j < numPot; j++) if (delR[j])
	{
		if (mPot == j) mPot = i;
		PotTB[i]->insertItem(PotTB[i]->count() - 1, PotTB[j]->currentText());
		PotTB[i]->setCurrentIndex(PotTB[i]->count() - 2);
		PotTB[j]->setCurrentIndex(PotTB[j]->count() - 1);
		if (Pot[j] != 0) 
		{
			disconnect(PotTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(PotTB[i], SIGNAL(DoubleClicked(int)), Pot[j], SLOT(show()));
		}
		else if (Pot[i] != 0)
		{
			disconnect(PotTB[i], SIGNAL(DoubleClicked(int)), Pot[i], SLOT(show()));
			connect(PotTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadPotential(int)));
		}
		pTable->setItem(i, 1, pTable->takeItem(j, 1));
		pTable->setItem(i, 2, pTable->takeItem(j, 2));
		Pot[i++] = Pot[j];
	}
	TBlock = false;
	numPot = i;
	refreshPotBox();
	if (mP != mPot) 
	{
		emit PotentialChanged();
		if (PotB != 0) PotB->setCurrentIndex(mPot);
	}
	if (Mol != 0) Mol->Changed();
	checkPotentials();
}

void ElState::deleteTermTables(bool *delR)
{
	int i, j, mT = mTerm;
	//printf("Beginn DeleteTermTable\n");
	TBlock = true;
	for (i=0; i < numTerm; i++) if (!delR[i])
	{
		//printf("deleteTermTable(%d)\n", i);
		if (mTerm == i) mTerm = 0;
		if (termTable[i] == 0) 
			disconnect(TermTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadTermTable(int)));
		else
		{
			disconnect(termTable[i], 0, this, SLOT(checkTermTables()));
			disconnect(TermTB[i], SIGNAL(DoubleClicked(int)), termTable[i], SLOT(show()));
			termTable[i]->setMolecule(0, 0);
		}
		termTable[i] = 0;
		if (tTable->item(i, 1) != 0) tTable->item(i, 1)->setText("");
		if (tTable->item(i, 2) != 0) tTable->item(i, 2)->setText("");
	}
	for (i=0; (i < numTerm ? delR[i] : false); i++) ;
	//printf("Vor Schleife\n");
	for (j=i+1; j < numTerm; j++) if (delR[j])
	{
		if (mTerm == j) mTerm = i;
		TermTB[i]->insertItem(TermTB[i]->count() - 1, TermTB[j]->currentText());
		TermTB[i]->setCurrentIndex(TermTB[i]->count() - 2);
		TermTB[j]->setCurrentIndex(TermTB[j]->count() - 1);
		if (termTable[j] != 0) 
		{
			disconnect(TermTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(TermTB[i], SIGNAL(DoubleClicked(int)), termTable[j], SLOT(show()));
		}
		else if (termTable[i] != 0)
		{
			disconnect(TermTB[i], SIGNAL(DoubleClicked(int)), termTable[i], SLOT(show()));
			connect(TermTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadTermTable(int)));
		}
		tTable->setItem(i, 1, tTable->takeItem(j, 1));
		tTable->setItem(i, 2, tTable->takeItem(j, 2));
		termTable[i++] = termTable[j];
	}
	numTerm = i;
	TBlock = false;
	refreshTermBox();
	if (mT != mTerm) 
	{
		emit TermTableChanged();
		if (TermB != 0) TermB->setCurrentIndex(mTerm);
	}
	if (Mol != 0) Mol->Changed();
}

void ElState::DunNameChanged()
{
	int i;
	QString N;
	for (i=0; i < numDun; i++) if (DunK[i] != 0) if ((N = DunK[i]->getName()) != DunTB[i]->currentText())
	{
		DunTB[i]->addItem(N);
		DunTB[i]->setCurrentIndex(DunTB[i]->count() - 1);
	}
	refreshDunBox();
}

void ElState::dunTableChanged(int r, int c)
{
	//printf("dunTablechanged(r=%d, c=%d)\n", r, c);
	if (r >= numDun || TBlock) return;
	if (DunK[r] == 0)
	{
		if (MW == 0) return;
		DunTable *D = MW->getDunTable(dTable->item(r, 1)->text(), Mol);
		if (D == 0) return;
		addDunTable(D, r);
	}
	if (c==1) DunK[r]->setFileName(dTable->item(r, c)->text());
	else if (c==2) DunK[r]->setSource(dTable->item(r, c)->text());
}

void ElState::dunTableClicked(int r, int c)
{
	//printf("dunTableClicked(r=%d, c=%d), numDun=%d\n", r, c, numDun);
	if (r != numDun || c!=0) return;
	if (DunTB == 0) DunTB = new ComboBox*[10];
	if (DunK == 0)
	{
		int n;
		DunK = new DunTable*[10];
		for (n=0; n<10; n++) DunK[n] = 0;
	}
	dTable->setCellWidget(r, c, DunTB[numDun] = new ComboBox(0));
	DunTB[numDun]->setNumber(numDun);
	connect(DunTB[numDun++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkDunTables()));
	refreshDunBox();
}

void ElState::FitDataNameChanged()
{
	int i;
	QString N;
	for (i=0; i < numFD; i++) if (FitD[i] != 0) if ((N = FitD[i]->getName()) != FDTB[i]->currentText())
	{
		FDTB[i]->addItem(N);
		FDTB[i]->setCurrentIndex(FDTB[i]->count() - 1);
	}
	refreshFitDataBox();
}

void ElState::fitDataTableChanged(int row, int column)
{
	if (row >= numFD || TBlock) return;
	if (FitD[row] == 0)
	{
		if (MW == 0) return;
		FitData *FD = MW->getFitData(fDTable->item(row, 1)->text(), Mol);
		if (FD == 0) return;
		addFitData(FD, row);
	}
	if (column == 1) FitD[row]->setFileName(fDTable->item(row, column)->text());
	else if (column == 2) FitD[row]->setSource(fDTable->item(row, column)->text());
}

void ElState::fitDataTableClicked(int row, int column)
{
	if (row != numFD || column != 0) return;
	if (FDTB == 0) FDTB = new ComboBox*[10];
	if (FitD == 0)
	{
		int n;
		FitD = new FitData*[10];
		for (n=0; n<10; n++) FitD[n] = 0;
	}
	fDTable->setCellWidget(row, column, FDTB[numFD] = new ComboBox(0));
	FDTB[numFD]->setNumber(numFD);
	connect(FDTB[numFD++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFitDataSets()));
	refreshFitDataBox();
}

double ElState::getBe()
{
    if (0 >= m_Be)
    {
        DunTable* DunT = getDunTable();
        if (0 != DunT)
        {
            DunT->getBe(m_Be);
            Mol->Changed();
        }
    }
    if (0 >= m_Be)
    {
        Potential* Pot = getPotential();
        if (0 != Pot)
        {
            m_Be = Pot->guessBe();
            Mol->Changed();
        }
    }
    return m_Be;
}

ComboBox *ElState::getDunBox()
{
	return DunB;
}

DunTable *ElState::getDunTable()
{
	return getDunTable(mDun);
}

QString ElState::getDunTableFileName(int i)
{
	if (i==-1) i = mDun;
	if (i<0 || i >= numDun) return "";
	return dTable->item(i, 1)->text();
}

DunTable *ElState::getDunTable(int i)
{
	if (i >= numDun) return 0;
	if (DunK[i] == 0 && MW != 0)
	{
		DunTable *D = MW->getDunTable(dTable->item(i, 1)->text(), Mol);
		if (D != 0) addDunTable(D, i);
	}
	return DunK[i];
}

QString ElState::getDunTableName(int i)
{
	if (i==-1) i = mDun;
	if (i<0 || i >= numDun) return "";
	return DunTB[i]->currentText();
}

QString ElState::getDunTableSource(int i)
{
	if (i==-1) i = mDun;
	if (i<0 || i >= mDun) return "";
	return dTable->item(i, 2)->text();
}

FitData* ElState::getFitData()
{
	return getFitData(mFD);
}

FitData* ElState::getFitData(int i)
{
	if (i >= numFD) return 0;
	if (FitD[i] == 0 && MW != 0) 
	{
		FitData *F = MW->getFitData(fDTable->item(i, 1)->text(), Mol);
		if (F != 0) addFitData(F, i);
	}
	return FitD[i];
}

ComboBox* ElState::getFitDataBox()
{
	return FDB;
}

QString ElState::getFitDataFileName(int i)
{
	if (i==-1) i = mFD;
	if (i<0 || i >= numFD) return "";
	return fDTable->item(i, 1)->text();
}

QString ElState::getFitDataName(int i)
{
	if (i==-1) i = mFD;
	if (i<0 || i >= numFD) return "";
	return FDTB[i]->currentText();
}

QString ElState::getFitDataSource(int i)
{
	if (i==-1) i = mFD;
	if (i<0 || i >= mFD) return "";
	return fDTable->item(i, 2)->text();
}

int ElState::getJStart(int Iso, int Comp) const
{
	int r = lambda;
	if (Mol->getJStep(Iso) == 2 && ((r == 0 && Par + Sym == 3) || (r == 1 && Par + Comp == 2))) r++;
	return r;
}

int ElState::getLambda()
{
	return lambda;
}

MainWindow *ElState::getMainWindow()
{
	return MW;
}

Molecule *ElState::getMolecule() const
{
	return Mol;
}

QString ElState::getName()
{
	return StateName->text();
}

int ElState::getNumDunTables()
{
	return numDun;
}

int ElState::getNumFitDataSets() const
{
	return numFD;
}

int ElState::getNumPotentials()
{
	return numPot;
}

int ElState::getNumTermTables()
{
	return numTerm;
}

int ElState::getOmega()
{
	return Omega;
}

int ElState::getParity()
{
	return Par;
}

ComboBox *ElState::getPotBox()
{
	return PotB;
}

Potential *ElState::getPotential()
{
	//printf("ElState::getPotential(): numPot=%d, mPot=%d\n", numPot, mPot);
	if (numPot == 0) return 0;
	if (Pot[mPot] == 0 && MW != 0) 
	{
		Potential *P = MW->getPotential(pTable->item(mPot, 1)->text(), Mol);
		if (P != 0) addPotential(P, mPot);
	}
	return Pot[mPot];
}

QString ElState::getPotentialFileName(int i)
{
	if (i == -1) i = mPot;
	if (i < 0 || i >= numPot) return "";
	return pTable->item(i, 1)->text();
}

Potential *ElState::getPotential(int i)
{
	if (Pot[i] == 0 && MW != 0)
	{
		Potential *P = MW->getPotential(pTable->item(i, 1)->text(), Mol);
		if (P != 0) addPotential(P, i);
	}
	return Pot[i];
}

QString ElState::getPotentialName(int i)
{
	if (i==-1) i = mPot;
	if (i < 0 || i >= numPot) return "";
	return PotTB[i]->currentText();
}

QString ElState::getPotentialSource(int i)
{
	if (i==-1) i = mPot;
	if (i<0 || i >= numPot) return "";
	return pTable->item(i, 2)->text();
}

float ElState::getS()
{
	return S = float(SB->currentIndex()) + (HS ? 0.5 : 0.0);
}

int ElState::getStateNum()
{
	return stateNum;
}

int ElState::getSymmetry()
{
	return Sym;
}

ComboBox *ElState::getTermBox()
{
	return TermB;
}

TermTable *ElState::getTermTable(bool calc)
{
	//printf("ElState::getTermTable(), sName=%s, mTerm=%d\n", getName().ascii(), mTerm);
	TermTable *T = 0;
	Potential *P;
	DunTable *D;
	if (numTerm != 0) 
	{
		if (termTable[mTerm] == 0 && MW != 0) 
			T = MW->getTermTable(tTable->item(mTerm, 1)->text(), Mol, this);
		if (T != 0) addTermTable(T, mTerm);
	}
	else if (MW == 0 || !calc) return 0;
	else if ((P = getPotential()) != 0) 
	{
		P->calcTermEnergies(T, 0, 0, false);
		addTermTable(T);
	}
	else if ((D = getDunTable()) != 0) 
	{
		D->calcTermEnergies(T, false);
		addTermTable(T);
	}
	else return 0;
	return termTable[mTerm];
}

QString ElState::getTermTableFileName(int i)
{
	if (i==-1) i = mTerm;
	if (i<0 || i >= numTerm) return "";
	return tTable->item(i, 1)->text();
}

TermTable *ElState::getTermTable(int i)
{
	if (termTable[i] == 0 && MW != 0)
	{
		TermTable *T = MW->getTermTable(tTable->item(i, 1)->text(), Mol, this);
		if (T != 0) addTermTable(T, i);
	}
	return termTable[i];
}

QString ElState::getTermTableName(int i)
{
	if (i==-1) i = mTerm;
	if (i<0 || i >= numTerm) return "";
	return TermTB[i]->currentText();
}

QString ElState::getTermTableSource(int i)
{
	if (i==-1) i = mTerm;
	if (i<0 || i >= numTerm) return "";
	return tTable->item(i, 2)->text();
}

void ElState::initialize(MainWindow *mw, Molecule *mol, int nS, int nlambda, int nOmega, int Parity,
						 int Symmetry, QString Letter)
{
	setMW(mw, mol);
	setS(nS);
	setLambda(nlambda);
	setOmega(nOmega);
	setParity(Parity);
	setSymmetry(Symmetry);
	setName(Letter);
}

bool ElState::isDunTableLoaded(int i)
{
	return (i >= 0 && i < numDun ? DunK[i] != 0 : (numDun > 0 ? DunK[mDun] != 0 : false));
}

bool ElState::isFitDataLoaded(int i)
{
	return (i >= 0 && i < numFD ? FitD[i] != 0 : (numFD > 0 ? FitD[mFD] != 0 : false));
}

bool ElState::isPotentialLoaded(int i)
{
	return (i >= 0 && i < numPot ? Pot[i] != 0 : (numPot > 0 ? Pot[mPot] != 0 : false));
}

bool ElState::isTermTableLoaded(int i)
{
	return (i >= 0 && i < numTerm ? termTable[i] != 0 
								  : (numTerm > 0 ? termTable[mTerm] != 0 : false));
}

void ElState::loadDunTable(int i)
{
	if (getDunTable(i) != 0) DunK[i]->show();
}

void ElState::loadFitData(int Number)
{
	if (getFitData(Number) != 0) FitD[Number]->show();
}

void ElState::loadPotential(int i)
{
	if (getPotential(i) != 0) Pot[i]->show();
}

void ElState::loadTermTable(int i)
{
	if (getTermTable(i) != 0) termTable[i]->show();
}

void ElState::PotNameChanged()
{
	int i;
	QString N;
	for (i=0; i < numPot; i++) if (Pot[i] != 0) if (PotTB[i]->currentText() != (N = Pot[i]->getName()))
	{
		PotTB[i]->addItem(N);
		PotTB[i]->setCurrentIndex(PotTB[i]->count() - 1);
	}
	refreshPotBox();
}

void ElState::potTableChanged(int r, int c)
{
	if (r >= numPot || TBlock) return;
	if (Pot[r] == 0)
	{
		if (MW == 0) return;
		Potential *P;
		if ((P = MW->getPotential(pTable->item(r, 1)->text())) == 0) return;
		addPotential(P, r);
	}
	if (c==1) Pot[r]->setFileName(pTable->item(r, c)->text());
	else if (c==2) Pot[r]->setSource(pTable->item(r, c)->text());
}

void ElState::potTableClicked(int r, int c)
{
	if (r != numPot || c!=0) return;
	if (PotTB == 0) PotTB = new ComboBox*[10];
	if (Pot == 0)
	{
		int n;
		Pot = new Potential*[10];
		for (n=0; n<10; n++) Pot[n] = 0;
	}
	pTable->setCellWidget(r, c, PotTB[numPot] = new ComboBox(0));
	PotTB[numPot]->setNumber(numPot);
	connect(PotTB[numPot++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkPotentials()));
	refreshPotBox();
}

void ElState::refreshDunBox()
{
	//printf("ElState::refreshDunBox(), numDun=%d\n", numDun);
	TBlock = true;
	int i, j, N = 0;
	QString T;
	DunTable *D;
	if (MW != 0) N = MW->getNumDunTables();
	//printf("S1\n");
	if (DunB != 0) 
	{
		DunB->clear();
		for (i=0; i < numDun; i++) DunB->addItem(DunTB[i]->currentText());
		DunB->setCurrentIndex(mDun);
	}
	//if (numDun > 0) printf("DunTB[0]=%d\n", DunTB[0]);
	for (i=0; i < numDun; i++)
	{
		//printf("i=%d, numDun=%d\n", i, numDun);
		T = DunTB[i]->currentText();
		//printf("T=%s\n", T.ascii());
		DunTB[i]->clear();
		//printf("nach clear\n");
		DunTB[i]->addItem(T);
		//printf("Ende Runde\n");
		DunTB[i]->addItem("");
	}
	//printf("S3\n");
	for (i=0; i<N; i++) if (!(D = MW->getDunTable(i))->isAssigned())
	{
		T = D->getName();
		for (j=0; j < numDun; j++) DunTB[j]->addItem(T);
		if (DunB != 0) DunB->addItem(T);
	}
	if (numDun==0 && DunB != 0) 
	{
		DunB->addItem("");
		DunB->setCurrentIndex(DunB->count() - 1);
	}
	TBlock = false;
	//printf("Ende\n");
}

void ElState::refreshFitDataBox()
{
	TBlock = true;
	int i, j, N = 0;
	QString T;
	FitData *FD;
	if (MW != 0) N = MW->getNumFitDataSets();
	if (FDB != 0)
	{
		FDB->clear();
		for (i=0; i < numFD; i++) FDB->addItem(FDTB[i]->currentText());
		FDB->setCurrentIndex(mFD);
	}
	for (i=0; i < numFD; i++)
	{
		T = FDTB[i]->currentText();
		FDTB[i]->clear();
		FDTB[i]->addItem(T);
		FDTB[i]->addItem("");
	}
	for (i=0; i<N; i++) if (!(FD = MW->getFitData(i))->isAssigned())
	{
		T = FD->getName();
		for (j=0; j < numFD; j++) FDTB[j]->addItem(T);
		if (FDB != 0) FDB->addItem(T);
	}
	if (numFD == 0 && FDB != 0)
	{
		FDB->addItem("");
		FDB->setCurrentIndex(FDB->count() - 1);
	}
	TBlock = false;
}

void ElState::refreshPotBox()
{
	TBlock = true;
	int i, j, N=0;
	QString T;
	Potential *P;
	if (MW != 0) N = MW->getNumPotentials();
	if (PotB != 0)
	{
		PotB->clear();
		for (i=0; i < numPot; i++) PotB->addItem(PotTB[i]->currentText());
		PotB->setCurrentIndex(mPot);
	}
	for (i=0; i < numPot; i++)
	{
		T = PotTB[i]->currentText();
		PotTB[i]->clear();
		PotTB[i]->addItem(T);
		PotTB[i]->addItem("");
	}
	for (i=0; i<N; i++) if (!(P = MW->getPotential(i))->isAssigned())
	{
		T = P->getName();
		for (j=0; j < numPot; j++) PotTB[j]->addItem(T);
		if (PotB != 0) PotB->addItem(T);
	}
	if (numPot == 0 && PotB != 0)
	{
		PotB->addItem("");
		PotB->setCurrentIndex(PotB->count() - 1);
	}
	TBlock = false;
}

void ElState::refreshTermBox()
{
	//printf("ElState::refreshTermBox, mTerm=%d, numTerm=%d\n", mTerm, numTerm);
	TBlock = true;
	int i, j, N=0;
	QString Text;
	TermTable *T;
	if (MW != 0) N = MW->getNumTermTables();
	if (TermB != 0)
	{
		TermB->clear();
		for (i=0; i < numTerm; i++) TermB->addItem(TermTB[i]->currentText());
		TermB->setCurrentIndex(mTerm);
	}
	for (i=0; i < numTerm; i++)
	{
		Text = TermTB[i]->currentText();
		TermTB[i]->clear();
		TermTB[i]->addItem(Text);
		TermTB[i]->addItem("");
	}
	for (i=0; i<N; i++) if (!(T = MW->getTermTable(i))->isAssigned())
	{
		Text = T->getName();
		for (j=0; j < numTerm; j++) TermTB[j]->addItem(Text);
		if (TermB != 0) TermB->addItem(Text);
	}
	if (numTerm == 0 && TermB != 0)
	{
		TermB->addItem("");
		TermB->setCurrentIndex(TermB->count() - 1);
	}
	TBlock = false;
}

void ElState::removeDunTable(int j)
{
	bool Del[numDun];
	int i;
	for (i=0; i < numDun; i++) Del[i] = true;
	Del[j] = false;
	deleteDunTables(Del);
}

void ElState::removeFitData(int j)
{
	bool Del[numFD];
	int i;
	for (i=0; i < numFD; i++) Del[i] = true;
	Del[j] = false;
	deleteFitDataSets(Del);
}

void ElState::removePotential(Potential* rPot)
{
	int i;
	for (i=0; i < numPot; i++) if (Pot[i] == rPot) removePotential(i);
}

void ElState::removePotential(int j)
{
	bool Del[numPot];
	int i;
	for (i=0; i < numPot; i++) Del[i] = true;
	Del[j] = false;
	deletePotentials(Del);
}

void ElState::removeTermTable(int j)
{
	bool Del[numTerm];
	int i;
	for (i=0; i < numTerm; i++) Del[i] = true;
	Del[j] = false;
	deleteTermTables(Del);
}

void ElState::setDunBox(ComboBox *Box)
{
	DunB = Box;
	connect(Box, SIGNAL(activated(QString)),
			this, SLOT(updateDunTable(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showDunTable()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showDunMenu(QPoint)));
}

void ElState::setFitDataBox(ComboBox* Box)
{
	FDB = Box;
	connect(Box, SIGNAL(activated(QString)), this, SLOT(updateFitData(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showFitData()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showFitDataMenu(QPoint)));
}

void ElState::setLambda(int nL)
{
	//printf("lambda=%d\n", nL);
	if (nL <= 0) lambdaB->setCurrentIndex(0);
	else if (nL == 1) lambdaB->setCurrentIndex(1);
	else lambdaB->setCurrentIndex(2);
}

void ElState::setMainDunTable(DunTable *D)
{
	mDun = numDun;
	addDunTable(D);
    if (!m_reading && 0 != D) D->getBe(m_Be);
}

void ElState::setMainDunTable(int i)
{
	mDun = i;
    if (!m_reading)
    {
        DunTable* D = getDunTable();
        if (0 != D) D->getBe(m_Be);
    }
	if (Mol != 0) Mol->Changed();
}

bool ElState::setMainDunTable(QString FileName)
{
	int i;
	for (i=0; (i < numDun ? dTable->item(i, 1)->text() != FileName : false); i++) ;
	if (i == numDun) return false;
    mDun = i;
    if (!m_reading)
    {
        DunTable* D = getDunTable();
        if (0 != D) D->getBe(m_Be);
    }
	return true;
}

void ElState::setMainFitData(FitData* set)
{
	mFD = numFD;
	addFitData(set);
}

bool ElState::setMainFitData(QString FileName)
{
	int i;
	for (i=0; (i < numFD ? fDTable->item(i, 1)->text() != FileName : false); i++) ;
	if (i == numFD) return false;
	mFD = i;
	return true;
}

void ElState::setMainFitData(int BoxIndex)
{
	mFD = BoxIndex;
	if (Mol != 0) Mol->Changed();
}

void ElState::setMainPotential(int i)
{
	mPot = i;
    if ((0 >= m_Be || 0 == numDun) && !m_reading)
    {
        Potential* Pot = getPotential();
        if (0 != Pot) m_Be = Pot->guessBe();
    }
	if (Mol != 0) Mol->Changed();
}

void ElState::setMainPotential(Potential *pot)
{
	mPot = numPot;
	addPotential(pot);
    if ((0 >= m_Be || 0 == numDun) && 0 != pot && !m_reading) m_Be = pot->guessBe();
}

bool ElState::setMainPotential(QString FileName)
{
	int i;
	for (i=0; (i < numPot ? pTable->item(i, 1)->text() != FileName : false); i++) ;
	//printf("ElState::setMainPotential(Filename): numPot=%d, i=%d\n", numPot, i);
	if (i == numPot) return false;
	mPot = i;
    if ((0 >= m_Be || 0 == numDun) && !m_reading)
    {
        Potential* Pot = getPotential();
        if (0 != Pot) m_Be = Pot->guessBe();
    }
	return true;
}

void ElState::setMainTermTable(int i)
{
	mTerm = i;
	if (Mol != 0) Mol->Changed();
}

bool ElState::setMainTermTable(QString FileName)
{
	int i;
	for (i=0; (i < numTerm ? tTable->item(i, 1)->text() != FileName : false); i++) ;
	if (i == numTerm) return false;
	mTerm = i;
	return true;
}

void ElState::setMainTermTable(TermTable *T)
{
	mTerm = numTerm;
	addTermTable(T);
}

void ElState::setMW(MainWindow *mw, Molecule *mol)
{
	//printf("ElState::setMW\n");
	MW = mw;
	Mol = mol;
	if (mol == 0) return;
	if (mol->getAtom1() != mol->getAtom2())
	{
		ParB->setCurrentIndex(0);
		ParB->setEnabled(false);
	}
	else ParB->setEnabled(true);
    if (mol->getAtom1() != 0 && mol->getAtom2() != 0 && 2 * ((mol->getAtom1()->getNP() + mol->getAtom2()->getNP()) / 2)
		!= mol->getAtom1()->getNP() + mol->getAtom2()->getNP())
	{
		SB->addItems(QStringList() << "1/2" << "3/2" << "5/2");
		HS = true;
	}
	else SB->addItems(QStringList() << "0" << "1" << "2");
	ShowDunAction = new QAction("Show Dunham table", Mol);
	ShowDunAction->setStatusTip("Show this table of Dunham coefficients");
	connect(ShowDunAction, SIGNAL(triggered()), this, SLOT(showDunTable()));
	ShowDMAction = new QAction("Show state properties", Mol);
	ShowDMAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowDMAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowDunMenu = new QMenu(Mol);
	ShowDunMenu->addAction(ShowDunAction);
	ShowDunMenu->addAction(ShowDMAction);
	ShowFitDataAction = new QAction("Show fit dataset", Mol);
	ShowFitDataAction->setStatusTip("Show this table of fit data");
	connect(ShowFitDataAction, SIGNAL(triggered()), this, SLOT(showFitData()));
	ShowFDMAction = new QAction("Show state properties", Mol);
	ShowFDMAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowFDMAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowFitDataMenu = new QMenu(Mol);
	ShowFitDataMenu->addAction(ShowFitDataAction);
	ShowFitDataMenu->addAction(ShowFDMAction);
	ShowPotAction = new QAction("Show potential", Mol);
	ShowPotAction->setStatusTip("Show this table of potential coeffecients/points");
	connect(ShowPotAction, SIGNAL(triggered()), this, SLOT(showPotential()));
	ShowPMAction = new QAction("Show state properties", Mol);
	ShowPMAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowPMAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowPotMenu = new QMenu(Mol);
	ShowPotMenu->addAction(ShowPotAction);
	ShowPotMenu->addAction(ShowPMAction);
	ShowTermAction = new QAction("Show term table", Mol);
	ShowTermAction->setStatusTip("Show this table of term energies");
	connect(ShowTermAction, SIGNAL(triggered()), this, SLOT(showTermTable()));
	ShowTMAction = new QAction("Show state properties", Mol);
	ShowTMAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowTMAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowTermMenu = new QMenu(Mol);
	ShowTermMenu->addAction(ShowTermAction);
	ShowTermMenu->addAction(ShowTMAction);
	ShowAction = new QAction("Show state properties", Mol);
	ShowAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowMenu = new QMenu(Mol);
	ShowMenu->addAction(ShowAction);
	connect(MW, SIGNAL(TermTablesChanged()), this, SLOT(refreshTermBox()));
	connect(MW, SIGNAL(DunTablesChanged()), this, SLOT(refreshDunBox()));
	connect(MW, SIGNAL(PotentialsChanged()), this, SLOT(refreshPotBox()));
	connect(MW, SIGNAL(FitDataChanged()), this, SLOT(refreshFitDataBox()));
	//printf("Ende SetMW\n");
}

void ElState::setName(QString N)
{
	StateName->setText(N);
	MDIChild::setName(N);
    if (Mol != 0)
    {
        Mol->Changed();
        Mol->setStateName(stateNum, N);
        setWindowTitle("State properties " + Mol->getName() + " " + N + ":");
    }
}

void ElState::setOmega(int O)
{
	int St = lambda - S;
	if (St < 0) St = 0;
	if (O >= St && O <= lambda + S) OmegaB->setCurrentIndex(O - St);
}

void ElState::setParity(int P)
{
	ParB->setCurrentIndex(P);
}

void ElState::setPotBox(ComboBox *Box)
{
	PotB = Box;
	connect(Box, SIGNAL(activated(QString)), this, SLOT(updatePotential(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showPotential()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showPotMenu(QPoint)));
}

void ElState::setS(float s)
{
	int i = int(S=s);
	if (i >= 0 && i < SB->count()) SB->setCurrentIndex(i);
}

void ElState::setStateNum(int N)
{
	stateNum = N;
}

void ElState::setSymmetry(int S)
{
	if (S == -1) S=2;
	if (S >= 0 && S <= 2) SymB->setCurrentIndex(S);
}

void ElState::setTermBox(ComboBox *Box)
{
	TermB = Box;
	connect(Box, SIGNAL(activated(QString)), 
			this, SLOT(updateTermTable(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showTermTable()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showTermMenu(QPoint)));
}

void ElState::showDunMenu(QPoint P)
{
	sDun = mDun;
	if (ShowDunMenu != 0) ShowDunMenu->popup(P);
}

void ElState::showDunTable()
{
	if (sDun >= 0 && sDun < numDun)
	{
		if (DunK[sDun] == 0) loadDunTable(sDun);
		else DunK[sDun]->show();
	}
}

void ElState::showFitData()
{
	if (sFD >= 0 && sFD < numFD)
	{
		if (FitD[sFD] == 0) loadFitData(sFD);
		else FitD[sFD]->show();
	}
}

void ElState::showFitDataMenu(QPoint P)
{
	sFD = mFD;
	if (ShowFitDataMenu != 0) ShowFitDataMenu->popup(P);
}

void ElState::showPotential()
{
	if (sPot >= 0 && sPot < numPot)
	{
		if (Pot[sPot] == 0) loadPotential(sPot);
		else Pot[sPot]->show();
	}
}

void ElState::showPotMenu(QPoint P)
{
	sPot = mPot;
	if (ShowPotMenu != 0) ShowPotMenu->popup(P);
}

void ElState::showShowMenu(QPoint P)
{
	if (ShowMenu != 0) ShowMenu->popup(P);
}

void ElState::showTermMenu(QPoint P)
{
	sTerm = mTerm;
	if (ShowTermMenu != 0) ShowTermMenu->popup(P);
}

void ElState::showTermTable()
{
	if (sTerm >= 0 && sTerm < numTerm)
	{
		if (termTable[sTerm] == 0) loadTermTable(sTerm);
		else termTable[sTerm]->show();
	}
}

void ElState::TermNameChanged()
{
	int i;
	QString N;
	for (i=0; i < numTerm; i++) if (termTable[i] != 0) 
			if ((N = termTable[i]->getName()) != TermTB[i]->currentText())
	{
		TermTB[i]->addItem(N);
		TermTB[i]->setCurrentIndex(TermTB[i]->count() - 1);
	}
	refreshTermBox();
}

void ElState::termTableChanged(int r, int c)
{
	if (r >= numTerm || TBlock) return;
	if (termTable[r] == 0)
	{
		if (MW == 0) return;
		TermTable *T;
		if ((T = MW->getTermTable(tTable->item(r, 1)->text(), Mol, this)) == 0) return;
		addTermTable(T, r);
	}
	if (c==1) termTable[r]->setFileName(tTable->item(r, c)->text());
	else if (c==2) termTable[r]->setSource(tTable->item(r, c)->text());
}

void ElState::termTableClicked(int r, int c)
{
	if (r != numTerm || c!=0) return;
	if (TermTB == 0) TermTB = new ComboBox*[10];
	if (termTable == 0) 
	{
		int n;
		termTable = new TermTable*[10];
		for (n=0; n<10; n++) termTable[n] = 0;
	}
	tTable->setCellWidget(r, c, TermTB[numTerm] = new ComboBox(0));
	TermTB[numTerm]->setNumber(numTerm);
	connect(TermTB[numTerm++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkTermTables()));
	refreshTermBox();
}

void ElState::TypeChanged()
{
	bool Changed = false;
	int St, i;
	if (lambdaB->currentIndex() != lambda)
	{
		Changed = true;
		lambda = lambdaB->currentIndex();
		if (lambda == 0) SymB->setEnabled(true);
		else
		{
			SymB->setCurrentIndex(0);
			SymB->setEnabled(false);
		}
		S = -1;
	}
	if (Omega != OmegaB->currentIndex())
	{
		Changed = true;
		St = lambda - S;
		if (St < 0) St = 0;
		Omega = OmegaB->currentIndex() + St;
	}
	if (S != SB->currentIndex())
	{
		Changed = true;
		S = SB->currentIndex();
		if (HS) S = 0.5 * (S+1);
		OmegaB->blockSignals(true);
		OmegaB->clear();
		for (i = St = lambda - S; i <= lambda + S; i++) OmegaB->addItem(QString::number(i));
		if (Omega < St) Omega = St;
		if (Omega > lambda + S) Omega = lambda + S;
		OmegaB->setCurrentIndex(Omega - St);
		if (S == 0) OmegaB->setEnabled(false);
		else OmegaB->setEnabled(true);
		OmegaB->blockSignals(false);
	}
	if (Par != ParB->currentIndex())
	{
		Changed = true;
		Par = ParB->currentIndex();
	}
	if (SymB->currentIndex() != (Sym == -1 ? 2 : Sym))
	{
		Changed = true;
		Sym = SymB->currentIndex();
		if (Sym == 2) Sym = -1;
	}
	if (StateName->text() != (Name = getName()))
	{
		Changed = true;
		MDIChild::setName(StateName->text());
		setWindowTitle("State properties " + Mol->getName() + " " + StateName->text() + ":");
		Mol->setStateName(stateNum, Name.left(1));
	}
	if (Changed) emit PropertiesChanged(stateNum);
	if (Mol != 0) Mol->Changed();
}

void ElState::updateDunTable(QString Name)
{
	if (MW == 0)
	{
		printf("ElState::updateDunTable: error MW = 0!\n");
		return;
	}
	int i = DunB->currentIndex();
	if (i < numDun) setMainDunTable(i);
	else
	{
		int N = MW->getNumDunTables();
		DunTable *D = 0;
		for (i=0; (i<N ? (D = MW->getDunTable(i))->getName() != Name : false); i++) ;
		if (i<N) setMainDunTable(D);
	}
	emit DunhamChanged();
}

void ElState::updateFitData(QString Name)
{
	if (MW == 0)
	{
		printf("ElState::updateFitData: error MW = 0!\n");
		return;
	}
	int i = FDB->currentIndex();
	if (i < numFD) setMainFitData(i);
	else
	{
		int N = MW->getNumFitDataSets();
		FitData *FD = 0;
		for (i=0; (i<N ? (FD = MW->getFitData(i))->getName() != Name : false); i++) ;
		if (i<N) setMainFitData(FD);
	}
	emit FitDataChanged();
}

void ElState::updatePotential(QString Name)
{
	//printf("ElState::updatePotential\n");
	if (MW == 0)
	{
		printf("ElState::updatePotential: error MW = 0!\n");
		return;
	}
	int i = PotB->currentIndex();
	if (i < numPot) setMainPotential(i);
	else
	{
		int N = MW->getNumPotentials();
		Potential *P = 0;
		for (i=0; (i<N ? (P = MW->getPotential(i))->getName() != Name : false); i++) ;
		if (i<N) setMainPotential(P);
	}
	//printf("Ende updatePotential\n");
	emit PotentialChanged();
}

void ElState::updateTermTable(QString Name)
{
	if (MW == 0) 
	{
		printf("ElState::updateTermTable: error MW = 0!\n");
		return;
	}
	int i = TermB->currentIndex();
	if (i < numTerm) setMainTermTable(i);
	else
	{
		int N = MW->getNumTermTables();
		TermTable *T = 0;
		for (i=0; (i<N? (T = MW->getTermTable(i))->getName() != Name : false); i++) ;
		if (i<N) setMainTermTable(T);
	}
	emit TermTableChanged();
}

bool ElState::writeData(QString FileName)
{
    if (Mol == 0) return MDIChild::writeData(FileName);
	if (QMessageBox::information(this, "MolSpektAnalysis", "The electronic state is part of the molecule " + Mol->getName() 
			+ " and connot be saved separately. Instead the data of the molecule will be written to disk.", 
			QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok) == QMessageBox::Ok)
		return Mol->writeData(FileName);
	return true;
}
