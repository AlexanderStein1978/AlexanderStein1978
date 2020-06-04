//
// C++ Implementation: fcftable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "fcftable.h"
#include "molecule.h"
#include "elstate.h"
#include "potential.h"
#include "constants.h"
#include "utils.h"
#include "termtable.h"

#include <stdio.h>

#include <QStringList>
#include <QBrush>
#include <QColor>
#include <QFont>

FCFTable::FCFTable(MainWindow *mw, Molecule *Mol, QString LState, QString UState) : TableWindow(MDIChild::FranckCondonView, mw)
{
	printf("FCFTable::FCFTable\n");
	int i;
	JDiff = 0;
	molecule = Mol;
	View->addItem("Franck-Condon factors");
	View->addItem("FCF * therm. population");
	TUnit->addItem(QString("°C").right(2));
	TUnit->addItem("K");
	Js->setText("0");
	Jss->setText("0");
	setWindowTitle("FCF table");
	if (mw != 0)
	{
		int i, N = mw->getNumMolecules();
		for (i=0; i<N; i++) 
		{
			MolBox->addItem((Mol = mw->getMolecule(i))->getName());
			if (Mol == molecule) MolBox->setCurrentIndex(i);
		}
	}
	UpdateMolData();
	i = USBox->findText(UState);
	if (i!=-1) USBox->setCurrentIndex(i);
	i = LSBox->findText(LState);
	if (i!=-1) LSBox->setCurrentIndex(i);
	connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateMolData()));
	connect(Js, SIGNAL(editingFinished()), this, SLOT(UpdateJ()));
	connect(Calc, SIGNAL(clicked()), this, SLOT(FillTable()));
	connect(View, SIGNAL(currentIndexChanged(int)), this, SLOT(ViewChanged(int)));
}

FCFTable::~FCFTable()
{
}

void FCFTable::FillTable()
{
	//printf("FCFTable::FillTable\n");
	char FCFF = 'f';
    int i, j, k, nvs = MaxvFCF, nvss = MaxvFCF, NumWFPoints = NumFCF_WFPoints;
	double ***FCF;
	QStringList HL;
	QString vBuff;
	QTableWidgetItem *Item;
	QFont FCFFont;
	TermTable *TT;
	FCFFont.setBold(true);
	if (molecule == 0) 
	{
		printf("FCFTable::FillTable: error molecule=0!\n");
		return;
	}
	Potential *uPot = molecule->getPot(USBox->currentText());
	Potential *lPot = molecule->getPot(LSBox->currentText());
	if (uPot == 0 || lPot == 0) return;
	lPot->getFCF(uPot, Iso->currentIndex(), i = Js->text().toInt(), j = Jss->text().toInt(), nvs, nvss,
                  FCF, NumWFPoints);
	JDiff = j - i;
	if (View->currentIndex() == 1)
	{
		double **Pop, T = Temp->text().toInt();
		int Nv, NJ;
		if (TUnit->currentIndex() == 0) T += C_K_C;
		ElState *State = lPot->getElState();
		if (State != 0) TT = State->getTermTable();
		if (TT != 0 ? TT->getMaxJ() < j || TT->getMaxv() < nvss - 1: true) 
			lPot->calcTermEnergies(TT, nvss, 10/*cMaxJ + 1*/, false);
			//printf("TT->maxJ=%d, TT->maxv=%d\n", TT->getMaxJ(), TT->getMaxv());
		TT->getThermPopulation(Pop, Nv, NJ, T, Iso->currentIndex(), 0);
		//printf("nvss=%d, j=%d, Nv=%d, NJ=%d\n", nvss, j, Nv, NJ);
		for (i=0; i < nvss; i++) for (k=0; k < nvs; k++) FCF[k][i][0] *= Pop[i][j];
		Destroy(Pop, Nv);
		FCFF = 'g';
	}
	Tab->setRowCount(2 * nvs);
	Tab->setColumnCount(nvss);
	//printf("nvss=%d, nvs=%d\n", nvss, nvs);
	for (i=0; i < nvss; i++) HL << "v''=" + QString::number(i);
	Tab->setHorizontalHeaderLabels(HL);
	HL.clear();
	for (j=0; j < nvs; j++) 
	{
		vBuff = "v'=" + QString::number(j);
		HL <<  vBuff + " FCF" << vBuff + " freq.";
	}
	Tab->setVerticalHeaderLabels(HL);
	for (i=0; i < nvss; i++) for (j=0, k=-1; j < nvs; j++)
	{ 
		if ((Item = Tab->item(++k, i)) == 0) 
		{
			Tab->setItem(k, i, Item = new QTableWidgetItem(QString::number(FCF[j][i][0], FCFF, 6)));
			Item->setForeground(QBrush(QColor(0, 0, 255)));
			Item->setFont(FCFFont);
		}
		else Item->setText(QString::number(FCF[j][i][0], FCFF, 6));
		if ((Item = Tab->item(++k, i)) == 0)
			Tab->setItem(k, i, Item = new QTableWidgetItem(QString::number(FCF[j][i][1], 'f', 4)));
		else Item->setText(QString::number(FCF[j][i][1], 'f', 4));
	}
	//printf("Vor Destroy FCF\n");
	Destroy(FCF, nvs, nvss);
	//printf("Ende FillTable()\n");
}

void FCFTable::UpdateJ()
{
	printf("FCFTable::UpdateJ\n");
	Jss->setText(QString::number(Js->text().toInt() + JDiff));
}

void FCFTable::UpdateMolData()
{
	//printf("FCFTable::UpdateMolData\n");
	molecule = MW->getMolecule(MolBox->currentText());
	Iso->clear();
	LSBox->clear();
	USBox->clear();
	if (molecule == 0) return;
	int i, N = molecule->getNumStates();
	QString Text;
	IsoTab *IT = molecule->getIso();
	for (i=0; i < IT->numIso; i++) 
		Iso->addItem(QString::number(IT->mNumIso1[i]) + *IT->chSymb1 + QString::number(IT->mNumIso2[i])
				+ *IT->chSymb2);
	Iso->setCurrentIndex(IT->refIso);
	//printf("NumStates=%d\n", N);
	for (i=0; i<N; i++) if (!molecule->getPotName(i).isEmpty())
	{
		LSBox->addItem(Text = molecule->getState(i));
		USBox->addItem(Text);
	}
	delete IT;
}

void FCFTable::ViewChanged(int i)
{
	if (i==0)
	{
		TUnit->setEnabled(false);
		Temp->setEnabled(false);
	}
	else
	{
		TUnit->setEnabled(true);
		Temp->setEnabled(true);
	}
}
