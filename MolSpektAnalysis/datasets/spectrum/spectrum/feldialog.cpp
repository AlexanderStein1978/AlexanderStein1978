//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "feldialog.h"
#include "MainWindow.h"
#include "molecule.h"
#include "termtable.h"
#include "isotab.h"

#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>


FELDialog::FELDialog(MainWindow* mw, double initMinH): QDialog(mw)
{
	int n, N = (MW = mw)->getNumMolecules();
	QGridLayout *L = new QGridLayout(this);
	setWindowTitle("Select settings for emission line search");
	L->addWidget(new QLabel("Molecule:", this), 0, 0);
	L->addWidget(MolBox = new QComboBox(this), 0, 1);
	L->addWidget(new QLabel("Excited state:", this), 1, 0);
	L->addWidget(StateBox = new QComboBox(this), 1, 1);
	L->addWidget(new QLabel("Isotopologue:", this), 2, 0);
	L->addWidget(IsoBox = new QComboBox(this), 2, 1);
	L->addWidget(new QLabel("Min line intensity:", this), 3, 0);
	L->addWidget(MinHEdit = new QLineEdit(QString::number(initMinH, 'g', 4), this), 3, 1);
	L->addWidget(new QLabel("Frequency tol:", this), 4, 0);
	L->addWidget(TolEdit = new QLineEdit("0.02", this), 4, 1);
	L->addWidget(new QLabel("Min lines per band:", this), 5, 0);
	L->addWidget(MLPBEdit = new QLineEdit("10", this), 5, 1);
	L->setRowMinimumHeight(6, 20);
	L->addWidget(OK = new QPushButton("OK", this), 7, 0);
	L->addWidget(Cancel = new QPushButton("Cancel", this), 7, 1);
	for (n=0; n<N; n++) if (DataA(Mol = MW->getMolecule(n))) MolBox->addItem(Mol->getName());
	MolBox->setEditable(false);
	StateBox->setEditable(false);
	IsoBox->setEditable(false);
	connect(MolBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(MolChanged(QString)));
	connect(StateBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(StateChanged(QString)));
	connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
	MolChanged(MolBox->currentText());
}

bool FELDialog::DataA(MainWindow* MW)
{
	int n, N = MW->getNumMolecules();
	for (n=0; n<N; n++) if (DataA(MW->getMolecule(n))) return true;
	return false;
}

bool FELDialog::DataA(Molecule* Mol)
{
	int n, N = Mol->getNumStates();
	for (n=1; n<N; n++) if (DataA(Mol, Mol->getStateP(n))) return true;
	return false;
}

bool FELDialog::DataA(Molecule* Mol, ElState* State)
{
	int n, N = Mol->getNumIso();
	for (n=0; n<N; n++) if (DataA(Mol, State, n)) return true;
	return false;
}

bool FELDialog::DataA(Molecule* Mol, ElState* State, int Iso)
{
	ElState *XState = Mol->getStateP(0);
	TermTable *Term = (XState != 0 ? XState->getTermTable() : 0);
	if (Term == 0) return false;
	int n, NIso = Term->getNumIso(), *IsoZ = Term->getIsoZ();
	for (n=0; (n < NIso ? IsoZ[n] != Iso: false); n++) ;
	if (n == NIso) return false;
	if ((Term = State->getTermTable()) == 0) return false;
	for (n=0, NIso = Term->getNumIso(), IsoZ = Term->getIsoZ(); IsoZ[n] != Iso; n++) ;
	return n != NIso;
}

int FELDialog::getIso()
{
	int n, i, I = IsoBox->currentIndex();
	for (n=i=0; true; n++) if (DataA(Mol, State, n)) 
	{
		if (i==I) return n;
		i++;
	}
}

double FELDialog::getMinH()
{
	return MinHEdit->text().toDouble();
}

int FELDialog::getMLPB()
{
	return MLPBEdit->text().toInt();
}

Molecule* FELDialog::getMol()
{
	return Mol;
}

ElState* FELDialog::getState()
{
	return State;
}

double FELDialog::getTol()
{
	return TolEdit->text().toDouble();
}

void FELDialog::MolChanged(QString Name)
{
	int n;
	for (n=0; (Mol = MW->getMolecule(n))->getName() != Name; n++) ;
	StateBox->blockSignals(true);
	StateBox->clear();
	for (n=1; n < Mol->getNumStates(); n++) if (DataA(Mol, Mol->getStateP(n))) StateBox->addItem(Mol->getState(n));
	StateBox->blockSignals(false);
	StateChanged(StateBox->currentText());
}

void FELDialog::StateChanged(QString Name)
{
	IsoTab *Iso = Mol->getIso();
	int n;
	for (n=1; Name != Mol->getState(n); n++) ;
	State = Mol->getStateP(n);
	IsoBox->blockSignals(true);
	IsoBox->clear();
	for (n=0; n < Mol->getNumIso(); n++) if (DataA(Mol, State, n)) IsoBox->addItem(Iso->getIsoName(n));
	IsoBox->blockSignals(false);
	delete Iso;
}
