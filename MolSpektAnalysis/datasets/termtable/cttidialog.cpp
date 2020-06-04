//
// C++ Implementation: CTTIDialog
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "cttidialog.h"


CTTIDialog::CTTIDialog(MainWindow* MW) : QDialog(MW)
{
	int n, N;
	mw = MW;
	mol = 0;
	QGridLayout *L = new QGridLayout(this);
	L->addWidget(new QLabel("Maximum v:", this), 0, 0);
	L->addWidget(Maxv = new QLineEdit(this), 0, 1);
	L->addWidget(new QLabel("Maximum J:", this), 1, 0);
	L->addWidget(MaxJ = new QLineEdit(this), 1, 1);
	L->addWidget(new QLabel("Molecule:", this), 2, 0);
	L->addWidget(MolBox = new QComboBox(this), 2, 1);
	L->setRowMinimumHeight(3, 20);
	L->addWidget(new QLabel("Electronic states:", this), 4, 0, 1, 2);
	L->addWidget(SList = new QListWidget(this), 5, 0, 3, 1);
	L->addWidget(SBox = new QComboBox(this), 5, 1);
	L->addWidget(Add = new QPushButton("Add", this), 6, 1);
	L->addWidget(Remove = new QPushButton("Remove", this), 7, 1);
	L->addWidget(MTables = new QRadioButton("Multiple tables", this), 8, 0);
	L->addWidget(STable = new QRadioButton("Single table", this), 8, 1);
	L->setRowMinimumHeight(9, 20);
	L->addWidget(OK = new QPushButton("OK", this), 10, 0);
	L->addWidget(Cancel = new QPushButton("Cancel", this), 10, 1);
	MTables->setChecked(true);
	connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(molChanged(int)));
	connect(Add, SIGNAL(clicked()), this, SLOT(add()));
	connect(Remove, SIGNAL(clicked()), this, SLOT(remove()));
	connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
	for (n=0, N = mw->getNumMolecules(); n<N; n++) 
		MolBox->addItem(mw->getMolecule(n)->getName());
}

void CTTIDialog::add()
{
	if (!SBox->currentText().isEmpty()) 
	{
		SList->addItem(SBox->currentText());
		SBox->removeItem(SBox->currentIndex());
	}
}

void CTTIDialog::getData(Molecule *&Mol, ElState **&States, int& NStates, int& Mv, int& MJ, bool &singleTab)
{
	if ((NStates = SList->count()) == 0 || mol == 0) return;
	int n, m, N = mol->getNumStates();
	States = new ElState*[NStates];
	for (n=0; n < NStates; n++) for (m=0; m<N; m++) 
		if (mol->getState(m) == SList->item(n)->text()) States[n] = mol->getStateP(m);
	Mv = Maxv->text().toInt();
	MJ = MaxJ->text().toInt();
	Mol = mol;
	singleTab = STable->isChecked();
}

void CTTIDialog::molChanged(int i)
{
	int n;
	mol = mw->getMolecule(i);
	SBox->clear();
	if (mol != 0) for (n=0; n < mol->getNumStates(); n++) SBox->addItem(mol->getState(n));
}

void CTTIDialog::remove()
{
	if (SList->currentIndex().isValid())
	{
		QListWidgetItem *I = SList->takeItem(SList->currentRow());
		SBox->addItem(I->text());
		delete I;
	}
}
