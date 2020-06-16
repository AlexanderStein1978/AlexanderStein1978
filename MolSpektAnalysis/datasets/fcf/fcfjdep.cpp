//
// C++ Implementation: FCFJDep
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2010 - 2016
//
// Copyright: See README file that comes with this source code
//
//


#include "fcfjdep.h"
#include "MainWindow.h"
#include "potential.h"
#include "DiagWindow.h"
#include "molecule.h"
#include "termtable.h"
#include "utils.h"
#include "isotab.h"

#include <stdio.h>

#include <QGridLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QComboBox>


FCFJDep::FCFJDep(MainWindow *MW)
{
	int n, N;
	setWindowTitle("Plot J dependent Franck-Condon factors");
	setMinimumSize(350, 200);
	setMaximumSize(350, 200);
	mw = MW;
	lPot = uPot = 0;
	lData = uData = 0;
	lDiag = 0;
	Mol = 0;
	aMJ = 350;
	amJ = 0;
	QGridLayout *L = new QGridLayout(this);
	QLabel *La = new QLabel("Molecule:", this);
	L->addWidget(La, 0, 0);
	MolB = new QComboBox(this);
	MolB->setEditable(false);
	if (MW != 0) 
		for(n=0, N = MW->getNumMolecules(); n<N; n++) MolB->addItem(MW->getMolecule(n)->getName()); 
	connect(MolB, SIGNAL(currentIndexChanged(int)), this, SLOT(molChanged(int)));
	L->addWidget(MolB, 0, 1);
	La = new QLabel("Isotopologue:", this);
	L->addWidget(La, 0, 2);
	IsoB = new QComboBox(this);
	IsoB->setEditable(false);
	connect(IsoB, SIGNAL(currentIndexChanged(int)), this, SLOT(lPotChanged()));
	connect(IsoB, SIGNAL(currentIndexChanged(int)), this, SLOT(uPotChanged()));
	L->addWidget(IsoB, 0, 3);
	La = new QLabel("Upper state:", this);
	L->addWidget(La, 1, 0);
	uPotB = new QComboBox(this);
	uPotB->setEditable(false);
	connect(uPotB, SIGNAL(currentIndexChanged(int)), this, SLOT(uPotChanged()));
	L->addWidget(uPotB, 1, 1);
	La = new QLabel("Lower state:", this);
	L->addWidget(La, 1, 2);
	lPotB = new QComboBox(this);
	lPotB->setEditable(false);
	connect(lPotB, SIGNAL(currentIndexChanged(int)), this, SLOT(lPotChanged()));
	L->addWidget(lPotB, 1, 3);
	QGridLayout *L2 = new QGridLayout;
	La = new QLabel("v':", this);
	L2->addWidget(La, 0, 0);
	vsB = new QComboBox(this);
	vsB->setEditable(false);
	connect(vsB, SIGNAL(currentIndexChanged(int)), this, SLOT(vChanged()));
	L2->addWidget(vsB, 0, 1);
	La = new QLabel("v'':", this);
	L2->addWidget(La, 0, 2);
	vssB = new QComboBox(this);
	vssB->setEditable(false);
	connect(vssB, SIGNAL(currentIndexChanged(int)), this, SLOT(vChanged()));
	L2->addWidget(vssB, 0, 3);
	L2->addWidget(new QLabel("J'':", this), 0, 4);
	JssB = new QComboBox(this);
	JssB->addItems(QStringList() << "J'-1" << "J'" << "J'+1");
	JssB->setCurrentIndex(1);
	JssB->setEditable(false);
	connect(JssB, SIGNAL(currentIndexChanged(int)), this, SLOT(vChanged()));
	L2->addWidget(JssB, 0, 5);
	L2->addWidget(new QLabel("Min J:", this), 1, 0);
	minJB = new QComboBox(this);
	minJB->setEditable(false);
	connect(minJB, SIGNAL(currentIndexChanged(int)), this, SLOT(minJChanged(int)));
	L2->addWidget(minJB, 1, 1);
	L2->addWidget(new QLabel("Max J:", this), 1, 2);
	maxJB = new QComboBox(this);
	maxJB->setEditable(false);
	connect(maxJB, SIGNAL(currentIndexChanged(int)), this, SLOT(maxJChanged(int)));
	L2->addWidget(maxJB, 1, 3);
	Fit = new QCheckBox("Fit function", this);
	Fit->setCheckState(Qt::Unchecked);
	L2->addWidget(Fit, 1, 4, 1, 2);
	L->addLayout(L2, 2, 0, 1, 4);
	L->setRowMinimumHeight(3, 20);
	newB = new QPushButton("New plot", this);
	connect(newB, SIGNAL(clicked(bool)), this, SLOT(newPlot()));
	L->addWidget(newB, 4, 0);
	addB = new QPushButton("Add to plot", this);
	connect(addB, SIGNAL(clicked(bool)), this, SLOT(addToPlot()));
	L->addWidget(addB, 4, 1);
	clearB = new QPushButton("Clear plot", this);
	connect(clearB, SIGNAL(clicked(bool)), this, SLOT(clearPlot()));
	L->addWidget(clearB, 4, 2);
	closeB = new QPushButton("Close", this);
	connect(closeB, SIGNAL(clicked(bool)), this, SLOT(close()));
	L->addWidget(closeB, 4, 3);
	molChanged(0);
}

void FCFJDep::addToPlot()
{
	if (mw != 0) 
	{
		DiagWindow *nD = mw->getActiveDiagWindow();
		if (nD != 0 ? nD->getType() == MDIChild::FCFJDependency : false) lDiag = nD;
	}
	if (lDiag == 0 && uPot == 0 && lPot == 0) return;
    int JD = JssB->currentIndex() - 1, vs = vsB->currentIndex(), vss = vssB->currentIndex(), NumWFPoints = NumPoints;
	int Jm = minJB->currentIndex() + (JD == -1 ? 1 : 0), I = IsoB->currentIndex(), r, J, NR;
	int JM = Jm + (NR = (maxJB->currentIndex() + 1)), Jss;
	double **Data = Create(NR, 2), lDat;
	for (r=0, J = Jm, Jss = Jm + JD; J < JM; J++, r++, Jss++)
	{
		Data[r][0] = J*(J+1);
		lDat = lData[I][vss][Jss];
		//printf("uData=%f, lData=%f, ", uData[I][vs][J], lData[I][vss][Jss]);
        lPot->getFastFCF(I, uPot, J, 0, uData[I][vs][J], 1, &lDat, &Jss, Data[r] + 1, NumWFPoints);
		//printf("Data[%d]=(%d;%f), uData=%f, lData=%f\n", 
			//   r, J*(J+1), Data[r][1], uData[I][vs][J], lData[I][vss][Jss]);
	}
	lDiag->addData(Data, NR);
	if (!lDiag->isVisible()) lDiag->show();
	lDiag->setFocus();
	Destroy(Data, NR);
}

void FCFJDep::clearPlot()
{
	if (mw != 0) 
	{
		DiagWindow *nD = mw->getActiveDiagWindow();
		if (nD != 0 ? nD->getType() == MDIChild::FCFJDependency : false) lDiag = nD;
	}
	if (lDiag == 0) return;
	lDiag->clearData();
	if (!lDiag->isVisible()) lDiag->show();
	lDiag->setFocus();
}

void FCFJDep::lPotChanged()
{
	int i, n=0, N = lPotB->currentIndex(), I = IsoB->currentIndex();
	ElState *S = 0;
	if (Mol != 0) for (i=0; n<=N && i < Mol->getNumStates(); i++)
	{
		S = Mol->getStateP(i);
		if (S->getTermTable() != 0 && S->getPotential() != 0) n++;
	}
	if (n<=N) S=0;
	if (S!=0)
	{
		TermTable *TT = S->getTermTable();
		double ****D = TT->getData();
		NIso = TT->getNumIso();
		Nvss = TT->getMaxv() + 1;
		NJss = TT->getMaxJ() + 1;
		lData = (D!=0 ? D[0] : 0);
	}
	lPot = (lData != 0 ? S->getPotential() : 0);
	if (lPot != 0 && I < NIso) for (n=0; (n < Nvss ? lData[I][n][2] != 0 : false); n++) ;
	else n=0;
	for (i = vssB->count(); i<n; i++) vssB->addItem(QString::number(i));
	for (i = vssB->count(); i>n; i--) vssB->removeItem(i-1);
}

void FCFJDep::minJChanged(int i)
{
	if (lData == 0 || uData == 0)
	{
		maxJB->clear();
		return;
	}
	int AMJ = aMJ, M, n, I = IsoB->currentIndex(), vs = vsB->currentIndex(), vss = vssB->currentIndex();
	int JD = JssB->currentIndex() - 1;
	if (JD == -1) i++;
	amJ = i;
	if (AMJ <= i) AMJ = i+1;
	for (M=i+1; (M < NJs && M + JD < NJss ? 
		lData[I][vss][M + JD] != 0.0 && uData[I][vs][M] != 0.0 : false); M++) ;
	if (AMJ >= M) AMJ = M-1;
	maxJB->clear();
	for (n=i+1; n<M; n++) maxJB->addItem(QString::number(n));
	maxJB->setCurrentIndex(AMJ - i - 1);
}

void FCFJDep::maxJChanged(int i)
{
	aMJ = amJ + i + 1;
}

void FCFJDep::molChanged(int i)
{
	lPotB->clear();
	uPotB->clear();
	IsoB->clear();
	if ((Mol = (mw != 0 ? mw->getMolecule(i) : 0)) == 0) return;
	IsoTab *Iso = Mol->getIso();
	ElState *S;
	QString SN;
	int n, N;
	for (n=0; n < Iso->numIso; n++) IsoB->addItem(Iso->getIsoName(n));
	IsoB->setCurrentIndex(Iso->refIso);
	delete Iso;
	for (n=0, N = Mol->getNumStates(); n<N; n++)
	{
		S = Mol->getStateP(n);
		if (S->getTermTable() != 0 && S->getPotential() != 0)
		{
			lPotB->addItem(SN = S->getName());
			uPotB->addItem(SN);
		}
	}
	if (uPotB->count() > 1) uPotB->setCurrentIndex(1);
}

void FCFJDep::newPlot()
{
	lDiag = (mw != 0 ? mw->CreateDiagWindow(MDIChild::FCFJDependency) : 0);
	if (lDiag == 0) return;
	lDiag->setUnits("J(J+1)", "FCF");
	addToPlot();
}

void FCFJDep::uPotChanged()
{
	int i, n=0, N = uPotB->currentIndex(), I = IsoB->currentIndex();
	ElState *S = 0;
	if (Mol != 0) for (i=0; n<=N && i < Mol->getNumStates(); i++)
	{
		S = Mol->getStateP(i);
		if (S->getTermTable() != 0 && S->getPotential() != 0) n++;
	}
	if (n<=N) S=0;
	if (S!=0)
	{
		TermTable *TT = S->getTermTable();
		double ****D = TT->getData();
		NIso = TT->getNumIso();
		Nvs = TT->getMaxv() + 1;
		NJs = TT->getMaxJ() + 1;
		uData = (D!=0 ? D[0] : 0);
	}
	uPot = (uData != 0 ? S->getPotential() : 0);
	if (uPot != 0 && I < NIso) for (n=0; (n < Nvs ? uData[I][n][2] != 0 : false); n++) ;
	else n=0;
	for (i = vsB->count(); i<n; i++) vsB->addItem(QString::number(i));
	for (i = vsB->count(); i>n; i--) vsB->removeItem(i-1);
}

void FCFJDep::vChanged()
{
	if (lData == 0 || uData == 0)
	{
		minJB->clear();
		maxJB->clear();
		return;
	}
	int AmJ = amJ, M, n, I = IsoB->currentIndex(), vs = vsB->currentIndex(), vss = vssB->currentIndex();
	int JD = JssB->currentIndex() - 1;
	for (M = (JD >= 0 ? 0 : 1); (M < NJs && M + JD < NJss ? 
		lData[I][vss][M + JD] != 0.0 && uData[I][vs][M] != 0.0 : false); M++) ;
	M--;
	if (AmJ >= M) AmJ = M-1;
	minJB->clear();
	for (n = (JD >= 0 ? 0 : 1); n<M; n++) minJB->addItem(QString::number(n));
	minJB->setCurrentIndex(JD >= 0 ? AmJ : AmJ - 1);
}
