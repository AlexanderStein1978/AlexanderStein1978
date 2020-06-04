//
// C++ Implementation: SelectLineTableDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "selectlinetabledialog.h"


SelectLineTableDialog::SelectLineTableDialog(MainWindow* mw): QDialog(mw)
{
	int n, N;
	MW = mw;
	Mol = 0;
	Tr = 0;
	QGridLayout *L = new QGridLayout(this);
	setWindowTitle("Select a linetable");
	L->addWidget(new QLabel("Molecule:", this), 0, 0);
	L->addWidget(MolB = new QComboBox(this), 0, 1);
	MolB->setEditable(false);
	for (n=0, N = MW->getNumMolecules(); n<N; n++) MolB->addItem(MW->getMolecule(n)->getName());
	L->addWidget(new QLabel("Transition:", this), 1, 0);
	L->addWidget(TranB = new QComboBox(this), 1, 1);
	TranB->setEditable(false);
	L->addWidget(new QLabel("Linetable:", this), 2, 0);
	L->addWidget(LTB = new QComboBox(this), 2, 1);
	LTB->setEditable(false);
	L->setRowMinimumHeight(3, 20);
	L->addWidget(OK = new QPushButton("OK", this), 4, 0);
	L->addWidget(Cancel = new QPushButton("Cancel", this), 4, 1);
	connect(MolB, SIGNAL(currentIndexChanged(int)), this, SLOT(MoleculeChanged(int)));
	connect(TranB, SIGNAL(currentIndexChanged(int)), this, SLOT(TransitionChanged(int)));
	connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
	if (N>0) MoleculeChanged(0);
}

LineTable* SelectLineTableDialog::getLineTable()
{
	int i;
	if (Tr != 0 && (i = LTB->currentIndex()) >= 0) return Tr->getLineTable(i);
	return 0;
}

void SelectLineTableDialog::MoleculeChanged(int i)
{
	int n, N=0;
	TranB->blockSignals(true);
	TranB->clear();
	if (i>=0)
	{
		Mol = MW->getMolecule(i);
		for (n=0, N = Mol->getNumTransitions(); n<N; n++) TranB->addItem(Mol->getTransitionP(n)->getName());
	}
	else Mol = 0;
	TranB->blockSignals(false);
	TransitionChanged(N>0 ? 0 : -1);
}

void SelectLineTableDialog::TransitionChanged(int i)
{
	int n, N;
	LTB->blockSignals(true);
	LTB->clear();
	if (i>=0)
	{
		Tr = Mol->getTransitionP(i);
		for (n=0, N = Tr->getNumLineTables(); n<N; n++) LTB->addItem(Tr->getLineTableName(n));
	}
	else Tr = 0;
	LTB->blockSignals(false);
}
