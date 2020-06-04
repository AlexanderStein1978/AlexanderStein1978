//
// C++ Implementation: FSPDiag
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "fspdiag.h"


FSPDiag::FSPDiag(MainWindow* mw, Progression Pr): QWidget(mw)
{
	int n, N = mw->getNumLineTables();
	P = new Progression;
	*P = Pr;
	MW = mw;
	setAttribute(Qt::WA_DeleteOnClose);
	setWindowTitle("Find similar progression");
	QGridLayout *L = new QGridLayout(this);
	QLabel *La = new QLabel("LineTable to search in:", this);
	L->addWidget(La, 0, 0, 1, 2);
	LTBox = new QComboBox(this);
	L->addWidget(LTBox, 1, 0, 1, 2);
	for (n=0; n<N; n++) LTBox->addItem(MW->getLineTable(n)->getName());
	L->setRowMinimumHeight(2, 20);
	SearchB = new QPushButton("Search", this);
	L->addWidget(SearchB, 3, 0);
	connect(SearchB, SIGNAL(clicked()), this, SLOT(search()));
	CloseB = new QPushButton("Close", this);
	L->addWidget(CloseB, 3, 1);
	connect(CloseB, SIGNAL(clicked()), this, SLOT(close()));
}

FSPDiag::~FSPDiag()
{
	if (P->L != 0) delete[] P->L;
	delete P;
}

void FSPDiag::search()
{
	MW->getLineTable(LTBox->currentIndex())->findSimilarProgression(*P);
}
