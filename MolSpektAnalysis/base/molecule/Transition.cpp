//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "Transition.h"
#include "molecule.h"
#include "fcftab.h"
#include "linetable.h"

#include <QMenu>
#include <QGridLayout>


Transition::Transition()
{
	lowerState = upperState = 0;
	Lines = 0;
	fcf = 0;
	Mol = 0;
	FCFB = LineB = USB = LSB = 0;
	FCFTB = LineTB = 0;
	MW = 0;
	ShowMenu = ShowUSMenu = ShowLSMenu = ShowLineMenu = ShowFCFMenu = 0;
	mLine = nLines = mFCF = nFCF = 0;
	transNum = 0;
	TBlock = false;
	QGridLayout *L = new QGridLayout(this);
	setWindowTitle("Transition");
	L->addWidget(new QLabel("Lower state:", this), 0, 0);
	L->addWidget(lSB = new ComboBox(this), 0, 1);
	lSB->setEditable(false);
	connect(lSB, SIGNAL(currentIndexChanged(int)), this, SLOT(lSBChanged(int)));
	L->addWidget(new QLabel("Upper state:", this), 0, 2);
	L->addWidget(uSB = new ComboBox(this), 0, 3);
	uSB->setEditable(false);
	connect(uSB, SIGNAL(currentIndexChanged(int)), this, SLOT(uSBChanged(int)));
	L->addWidget(new QLabel("Strength:", this), 1, 0);
	L->addWidget(transStrength = new QLineEdit("0", this), 1, 1);
	QDoubleValidator *Valid = new QDoubleValidator(transStrength);
	Valid->setBottom(0.0);
	transStrength->setValidator(Valid);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(new QLabel("Line tables:", this), 3, 0);
	L->addWidget(MergeTables = new QPushButton("Merge tables", this), 3, 3);
	connect(MergeTables, SIGNAL(clicked()), this, SLOT(mergeLineTables()));
	L->addWidget(lTable = new Table(nLineR = 10, 3, this), 4, 0, 1, 4);
	lTable->setHorizontalHeaderLabels(QStringList() << "Name" << "File name" << "Source");
	connect(lTable, SIGNAL(cellChanged(int, int)), this, SLOT(lineTableChanged(int, int)));
	connect(lTable, SIGNAL(RowDeleted(bool*)), this, SLOT(deleteLineTables(bool*)));
	connect(lTable, SIGNAL(cellClicked(int, int)), this, SLOT(lineTableClicked(int, int)));
	L->setRowMinimumHeight(5, 20);
	L->addWidget(new QLabel("FCF Tables:", this), 6, 0);
	L->addWidget(FCFT = new Table(nFCFR = 10, 3, this), 7, 0, 1, 4);
	FCFT->setHorizontalHeaderLabels(QStringList() << "Name" << "File name" << "Source");
	connect(FCFT, SIGNAL(cellChanged(int,int)), this, SLOT(fcfTableChanged(int,int)));
	connect(FCFT, SIGNAL(RowDeleted(bool*)), this, SLOT(deleteFCFTables(bool*)));
	connect(FCFT, SIGNAL(cellClicked(int,int)), this, SLOT(fcfTableClicked(int,int)));
    updateName();
}

Transition::~Transition()
{
	printf("Transition::~Transition\n");
	if (LineTB != 0) delete[] LineTB;
	if (Lines != 0) delete[] Lines;
	if (FCFTB != 0) delete[] FCFTB;
	if (fcf != 0) delete[] fcf;
	printf("Ende ~Transition\n");
}

void Transition::addFCFTable(FCFTab* Tab)
{
	addFCFTable(Tab, -1);
}

void Transition::addFCFTable(QString Name, QString FileName, QString Source)
{
	if (FCFB == 0 && Mol != 0) Mol->TTransitionsClicked(transNum, 2);
	TBlock = true;
	int n;
	if (fcf == 0) 
	{
		fcf = new FCFTab*[10];
		for (n=0; n < 10; n++) fcf[n] = 0;
		emit Added(transNum, 2);
	}
	if (FCFTB == 0) FCFTB = new ComboBox*[10];
	if (nFCF == nFCFR)
	{
		FCFTab **nfcf = new FCFTab*[nFCFR + 10];
		ComboBox **nFCFTB = new ComboBox*[nFCFR + 10];
		for (n=0; n < nFCFR; n++)
		{
			nfcf[n] = fcf[n];
			nFCFTB[n] = FCFTB[n];
		}
		while (n < nFCFR + 10) nfcf[n++] = 0;
		delete[] fcf;
		delete[] FCFTB;
		fcf = nfcf;
		FCFTB = nFCFTB;
		FCFT->setRowCount(nFCFR += 10);
	}
	FCFT->setCellWidget(nFCF, 0, FCFTB[nFCF] = new ComboBox(0));
	FCFTB[nFCF]->addItem(Name);
	connect(FCFTB[nFCF], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFCFTables()));
	connect(FCFTB[nFCF], SIGNAL(DoubleClicked(int)), this, SLOT(loadFCFTable(int)));
	FCFTB[nFCF]->setNumber(nFCF);
	FCFT->setItem(nFCF, 1, new QTableWidgetItem(FileName));
	FCFT->setItem(nFCF, 2, new QTableWidgetItem(Source));
	fcf[nFCF++] = 0;
	refreshFCFBox();
	TBlock = false;
}

void Transition::addFCFTable(FCFTab* Tab, int Index)
{
	QString cName, nName;
	int n;
	if (FCFB == 0 && Mol != 0) Mol->TTransitionsClicked(transNum, 2);
	TBlock = true;
	if (Index >= 0 && Index < nFCF) cName = FCFTB[Index]->currentText();
	if (fcf == 0) 
	{
		fcf = new FCFTab*[10];
		for (n=0; n < 10; n++) fcf[n] = 0;
		emit Added(transNum, 3);
	}
	if (FCFTB == 0) FCFTB = new ComboBox*[10];
	if (Index > nFCF) Index = nFCF;
	if (nFCF == nFCFR)
	{
		FCFTab **nfcf = new FCFTab*[nFCFR + 10];
		ComboBox **nFCFTB = new ComboBox*[nFCFR + 10];
		for (n=0; n  < nFCFR; n++)
		{
			nfcf[n] = fcf[n];
			nFCFTB[n] = FCFTB[n];
		}
		while (n < nFCF + 10) nfcf[n++] = 0;
		delete[] fcf;
		delete[] FCFTB;
		fcf = nfcf;
		FCFTB = nFCFTB;
		FCFT->setRowCount(nFCFR += 10);
	}
	if (Index < 0)
	{
		Index = nFCF;
		FCFT->setCellWidget(nFCF, 0, FCFTB[nFCF] = new ComboBox(0));
		if (Tab != 0)
		{
			FCFTB[nFCF]->addItem(Tab->getName());
			FCFTB[nFCF]->setCurrentIndex(FCFTB[nFCF]->count() - 1);
		}
		connect(FCFTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFCFTables()));
	}
	else if (fcf[Index] != 0)
	{
		disconnect(fcf[Index], 0, this, SLOT(checkFCFTables()));
		disconnect(FCFTB[Index], SIGNAL(DoubleClicked(int)), fcf[Index], SLOT(show()));
		fcf[Index]->setMolecule(0);
	}
	else disconnect(FCFTB[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadFCFTable(int)));
	fcf[Index] = Tab;
	if (Tab != 0)
	{
		FCFTB[Index]->addItem(nName = Tab->getName());
		FCFTB[Index]->setCurrentIndex(FCFTB[Index]->count() - 1);
		FCFT->setItem(Index, 1, new QTableWidgetItem(Tab->getFileName()));
		FCFT->setItem(Index, 2, new QTableWidgetItem(Tab->getSource()));
		Tab->setMolecule(Mol);
		Tab->setTransition(this);
		Tab->setStates(upperState, lowerState);
		connect(FCFTB[Index], SIGNAL(DoubleClicked(int)), Tab, SLOT(show()));
		connect(Tab, SIGNAL(fileNameChanged()), this, SLOT(checkFCFTables()));
		connect(Tab, SIGNAL(nameChanged(QString)), this, SLOT(FCFNameChanged()));
		connect(Tab, SIGNAL(SourceChanged()), this, SLOT(checkFCFTables()));
	}
	if (Index == nFCF) nFCF++;
	//printf("Vor refreshFCFBox\n");
	refreshFCFBox();
	//printf("Ende von Transition::addFCFTable\n");
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void Transition::addLineTable(LineTable *Tab, int Index)
{
	QString cName, nName;
	int n;
	if (LineB == 0 && Mol != 0) Mol->TTransitionsClicked(transNum, 2);
	TBlock = true;
	if (Index >= 0 && Index < nLines) cName = LineTB[Index]->currentText();
	if (Lines == 0) 
	{
		Lines = new LineTable*[10];
		for (n=0; n < nLineR; n++) Lines[n] = 0;
		emit Added(transNum, 2);
	}
	if (LineTB == 0) LineTB = new ComboBox*[10];
	if (Index > nLines) Index = nLines;
	if (nLines == nLineR)
	{
		LineTable **newLines = new LineTable*[nLineR + 10];
		ComboBox **nLineTB = new ComboBox*[nLineR + 10];
		for (n=0; n < nLineR; n++)
		{
			newLines[n] = Lines[n];
			nLineTB[n] = LineTB[n];
		}
		while (n < nLineR + 10) newLines[n++] = 0;
		delete[] Lines;
		delete[] LineTB;
		Lines = newLines;
		LineTB = nLineTB;
		lTable->setRowCount(nLineR += 10);
	}
	if (Index < 0)
	{
		Index = nLines;
		lTable->setCellWidget(nLines, 0, LineTB[nLines] = new ComboBox(0));
		if (Tab != 0)
		{
			LineTB[nLines]->addItem(Tab->getName());
			LineTB[nLines]->setCurrentIndex(LineTB[nLines]->count() - 1);
		}
		connect(LineTB[Index], SIGNAL(currentIndexChanged(int)), this, SLOT(checkLineTables()));
	}
	else if (Lines[Index] != 0)
	{
		disconnect(Lines[Index], 0, this, SLOT(checkLineTables()));
		disconnect(LineTB[Index], SIGNAL(DoubleClicked(int)), Lines[Index], SLOT(show()));
		Lines[Index]->setMolecule(0);
	}
	else disconnect(LineTB[Index], SIGNAL(DoubleClicked(int)), this, SLOT(loadLineTable(int)));
	Lines[Index] = Tab;
	if (Tab != 0)
	{
		LineTB[Index]->addItem(nName = Tab->getName());
		LineTB[Index]->setCurrentIndex(LineTB[Index]->count() - 1);
		lTable->setItem(Index, 1, new QTableWidgetItem(Tab->getFileName()));
		lTable->setItem(Index, 2, new QTableWidgetItem(Tab->getSource()));
		Tab->setMolecule(Mol);
		Tab->setTransition(this);
		connect(LineTB[Index], SIGNAL(DoubleClicked(int)), Tab, SLOT(show()));
		connect(Tab, SIGNAL(fileNameChanged()), this, SLOT(checkLineTables()));
		connect(Tab, SIGNAL(nameChanged(QString)), this, SLOT(LineNameChanged()));
		connect(Tab, SIGNAL(SourceChanged()), this, SLOT(checkLineTables()));
	}
	if (Index == nLines) nLines++;
	//printf("Vor refreshLineBox\n");
	refreshLineBox();
	//printf("Ende von Transition::addLineTable\n");
	TBlock = false;
	if (Mol != 0 && cName != nName) Mol->Changed();
}

void Transition::addLineTable(QString Name, QString FileName, QString Source)
{
	if (LineB == 0 && Mol != 0) Mol->TTransitionsClicked(transNum, 2);
	TBlock = true;
	int n;
	if (Lines == 0) 
	{
		Lines = new LineTable*[10];
		for (n=0; n < 10; n++) Lines[n] = 0;
		emit Added(transNum, 2);
	}
	if (LineTB == 0) LineTB = new ComboBox*[10];
	if (nLines == nLineR)
	{
		LineTable **newLines = new LineTable*[nLineR + 10];
		ComboBox **nLineTB = new ComboBox*[nLineR + 10];
		for (n=0; n < nLineR; n++)
		{
			newLines[n] = Lines[n];
			nLineTB[n] = LineTB[n];
		}
		while (n < nLineR + 10) newLines[n++] = 0;
		delete[] Lines;
		delete[] LineTB;
		Lines = newLines;
		LineTB = nLineTB;
		lTable->setRowCount(nLineR += 10);
	}
	lTable->setCellWidget(nLines, 0, LineTB[nLines] = new ComboBox(0));
	LineTB[nLines]->addItem(Name);
	connect(LineTB[nLines], SIGNAL(currentIndexChanged(int)), this, SLOT(checkLineTables()));
	connect(LineTB[nLines], SIGNAL(DoubleClicked(int)), this, SLOT(loadLineTable(int)));
	LineTB[nLines]->setNumber(nLines);
	lTable->setItem(nLines, 1, new QTableWidgetItem(FileName));
	lTable->setItem(nLines, 2, new QTableWidgetItem(Source));
	Lines[nLines++] = 0;
	refreshLineBox();
	TBlock = false;
}

void Transition::addLineTable(LineTable *Tab)
{
	addLineTable(Tab, -1);
}

void Transition::checkFCFTables()
{
	if (MW == 0 || TBlock) return;
	int i, d, D = MW->getNumFCFTables();
	bool Changed = false, TC = false;
	QString N, Text;
	FCFTab *FT = 0;
	QTableWidgetItem *I;
	for (i=0; i < nFCF; i++)
	{
		N = FCFTB[i]->currentText();
		if (N.isEmpty())
		{
			bool Del[nFCF];
			for (d=0; d < nFCF; d++) Del[d] = true;
			Del[i] = false;
			deleteFCFTables(Del);
			return;
		}
		else
		{
			if (fcf[i] != 0 ? N != fcf[i]->getName() : FCFTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (FT = MW->getFCFTable(d))->getName() != N : false); d++) ;
				if (d<D) 
				{
					addFCFTable(FT, i);
					Changed = true;
					TC = true;
				}
			}
			if (fcf[i] != 0)
			{
				if ((I = FCFT->item(i, 1))->text() != (Text = fcf[i]->getFileName())) 
				{
					I->setText(Text);
					TC = true;
				}
				if ((I = FCFT->item(i, 2))->text() != (Text = fcf[i]->getSource())) 
				{
					I->setText(Text);
					TC = true;
				}
			}
		}
	}
	if (Changed) refreshFCFBox();
	if (TC && Mol != 0) Mol->Changed();
}

void Transition::checkLineTables()
{
	if (MW == 0 || TBlock) return;
	int i, d, D = MW->getNumLineTables();
	bool Changed = false, TC = false;
	QString N, Text;
	LineTable *LT = 0;
	QTableWidgetItem *I;
	for (i=0; i < nLines; i++)
	{
		N = LineTB[i]->currentText();
		if (N.isEmpty())
		{
			bool Del[nLines];
			for (d=0; d < nLines; d++) Del[d] = true;
			Del[i] = false;
			deleteLineTables(Del);
			return;
		}
		else
		{
			if (Lines[i] != 0 ? N != Lines[i]->getName() : LineTB[i]->currentIndex() > 0)
			{
				for (d=0; (d<D ? (LT = MW->getLineTable(d))->getName() != N : false); d++) ;
				if (d<D) 
				{
					addLineTable(LT, i);
					Changed = true;
					TC = true;
				}
			}
			if (Lines[i] != 0)
			{
				if ((I = lTable->item(i, 1))->text() != (Text = Lines[i]->getFileName())) 
				{
					I->setText(Text);
					TC = true;
				}
				if ((I = lTable->item(i, 2))->text() != (Text = Lines[i]->getSource())) 
				{
					I->setText(Text);
					TC = true;
				}
			}
		}
	}
	if (Changed) refreshLineBox();
	if (TC && Mol != 0) Mol->Changed();
}

void Transition::deleteFCFTables(bool* delR)
{
	int i, j, mF = mFCF;
	TBlock = true;
	for (i=0; i < nFCF; i++) if (!delR[i])
	{
		if (mFCF == i) mFCF = 0;
		if (fcf[i] == 0) 
			disconnect(FCFTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadFCFTable(int)));
		else
		{
			disconnect(fcf[i], 0, this, SLOT(checkFCFTables()));
			disconnect(FCFTB[i], SIGNAL(DoubleClicked(int)), fcf[i], SLOT(show()));
			fcf[i]->setMolecule(0);
		}
		fcf[i] = 0;
		if (FCFT->item(i, 1) != 0) FCFT->item(i, 1)->setText("");
		if (FCFT->item(i, 2) != 0) FCFT->item(i, 2)->setText("");
		FCFTB[i]->setCurrentIndex(FCFTB[i]->count() - 1);
	}
	for (i=0; (i < nFCF ? delR[i] : false); i++) ;
	for (j=i+1; j < nFCF; j++) if (delR[j])
	{
		if (mFCF == j) mFCF = i;
		FCFTB[i]->insertItem(FCFTB[i]->count() - 1, FCFTB[j]->currentText());
		FCFTB[i]->setCurrentIndex(FCFTB[i]->count() - 2);
		FCFTB[j]->setCurrentIndex(FCFTB[j]->count() - 1);
		if (fcf[j] != 0)
		{
			disconnect(FCFTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(FCFTB[i], SIGNAL(DoubleClicked(int)), fcf[j], SLOT(show()));
		}
		else if (fcf[i] != 0)
		{
			disconnect(FCFTB[i], SIGNAL(DoubleClicked(int)), fcf[i], SLOT(show()));
			connect(FCFTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadFCFTable(int)));
		}
		FCFT->setItem(i, 1, FCFT->takeItem(j, 1));
		FCFT->setItem(i, 2, FCFT->takeItem(j, 2));
		fcf[i++] = fcf[j];
	}
	nFCF = i;
	TBlock = false;
	refreshFCFBox();
	if (mF != mFCF) 
	{
		emit FCFTableChanged();
		if (FCFB != 0) FCFB->setCurrentIndex(mFCF);
	}
	if (Mol != 0) Mol->Changed();
}

void Transition::deleteLineTables(bool *delR)
{
	int i, j, mL = mLine;
	TBlock = true;
	for (i=0; i < nLines; i++) if (!delR[i])
	{
		if (mLine == i) mLine = 0;
		if (Lines[i] == 0) 
			disconnect(LineTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadLineTable(int)));
		else
		{
			disconnect(Lines[i], 0, this, SLOT(checkLineTables()));
			disconnect(LineTB[i], SIGNAL(DoubleClicked(int)), Lines[i], SLOT(show()));
			Lines[i]->setMolecule(0);
		}
		Lines[i] = 0;
		if (lTable->item(i, 1) != 0) lTable->item(i, 1)->setText("");
		if (lTable->item(i, 2) != 0) lTable->item(i, 2)->setText("");
		LineTB[i]->setCurrentIndex(LineTB[i]->count() - 1);
	}
	for (i=0; (i < nLines ? delR[i] : false); i++) ;
	for (j=i+1; j < nLines; j++) if (delR[j])
	{
		if (mLine == j) mLine = i;
		LineTB[i]->insertItem(LineTB[i]->count() - 1, LineTB[j]->currentText());
		LineTB[i]->setCurrentIndex(LineTB[i]->count() - 2);
		LineTB[j]->setCurrentIndex(LineTB[j]->count() - 1);
		if (Lines[j] != 0)
		{
			disconnect(LineTB[i], SIGNAL(DoubleClicked(int)), 0, 0);
			connect(LineTB[i], SIGNAL(DoubleClicked(int)), Lines[j], SLOT(show()));
		}
		else if (Lines[i] != 0)
		{
			disconnect(LineTB[i], SIGNAL(DoubleClicked(int)), Lines[i], SLOT(show()));
			connect(LineTB[i], SIGNAL(DoubleClicked(int)), this, SLOT(loadLineTable(int)));
		}
		lTable->setItem(i, 1, lTable->takeItem(j, 1));
		lTable->setItem(i, 2, lTable->takeItem(j, 2));
		Lines[i++] = Lines[j];
	}
	nLines = i;
	TBlock = false;
	refreshLineBox();
	if (mL != mLine) 
	{
		emit LineTableChanged();
		if (LineB != 0) LineB->setCurrentIndex(mLine);
	}
	if (Mol != 0) Mol->Changed();
}

ComboBox *Transition::getFCFBox()
{
	return FCFB;
}

ComboBox *Transition::getLineBox()
{
	return LineB;
}

FCFTab *Transition::getFCFTable()
{
	return getFCFTable(mFCF);
}

LineTable *Transition::getLineTable()
{
	return getLineTable(mLine);
}

QString Transition::getFCFTableFileName(int i)
{
	if (i==-1) i = mFCF;
	if (i<0 || i >= nFCF) return "";
	return FCFT->item(i, 1)->text();
}

QString Transition::getLineTableFileName(int i)
{
	if (i==-1) i = mLine;
	if (i<0 || i >= nLines) return "";
	return lTable->item(i, 1)->text();
}

FCFTab *Transition::getFCFTable(int i)
{
	if (i >= nFCF) return 0;
	if (fcf[i] == 0 && MW != 0) 
	{
		FCFTab *F = MW->getFCFTable(FCFT->item(i, 1)->text(), Mol);
		if (F != 0) addFCFTable(F, i);
	}
	return fcf[i];
}

LineTable *Transition::getLineTable(int i)
{
	//printf("Transition::getLineTable: i=%d, mLine=%d\n", i, mLine);
	if (i >= nLines) return 0;
	if (Lines[i] == 0 && MW != 0) 
	{
		LineTable *L = MW->getLineTable(lTable->item(i, 1)->text(), Mol);
		if (L != 0) addLineTable(L, i);
	}
	return Lines[i];
}

QString Transition::getFCFTableName(int i)
{
	if (i==-1) i = mFCF;
	if (i<0 || i >= nFCF) return "";
	return FCFTB[i]->currentText();
}

QString Transition::getLineTableName(int i)
{
	if (i==-1) i = mLine;
	if (i<0 || i >= nLines) return "";
	return LineTB[i]->currentText();
}

QString Transition::getFCFTableSource(int i)
{
	if (i==-1) i = mFCF;
	if (i<0 || i >= nFCF) return "";
	return FCFT->item(i, 2)->text();
}

QString Transition::getLineTableSource(int i)
{
	if (i==-1) i = mLine;
	if (i<0 || i >= nLines) return "";
	return lTable->item(i, 2)->text();
}

ElState *Transition::getLowerState()
{
	return lowerState;
}

ComboBox *Transition::getLSB()
{
	return LSB;
}

MainWindow *Transition::getMainWindow()
{
	return MW;
}

Molecule *Transition::getMolecule()
{
	return Mol;
}

int Transition::getNumFCFTables()
{
	return nFCF;
}

int Transition::getNumLineTables()
{
	return nLines;
}

ElState *Transition::getUpperState()
{
	return upperState;
}

ComboBox *Transition::getUSB()
{
	return USB;
}

void Transition::FCFNameChanged()
{
	int i;
	QString N;
	for (i=0; i < nFCF; i++) if (fcf[i] != 0) if ((N = fcf[i]->getName()) != FCFTB[i]->currentText())
	{
		FCFTB[i]->addItem(N);
		FCFTB[i]->setCurrentIndex(FCFTB[i]->count() - 1);
	}
	refreshFCFBox();
}

void Transition::LineNameChanged()
{
	int i;
	QString N;
	for (i=0; i < nLines; i++) if (Lines[i] != 0) if ((N = Lines[i]->getName()) != LineTB[i]->currentText())
	{
		LineTB[i]->addItem(N);
		LineTB[i]->setCurrentIndex(LineTB[i]->count() - 1);
	}
	refreshLineBox();
}

void Transition::fcfTableChanged(int r, int c)
{
	if (r >= nFCF || TBlock) return;
	if (fcf[r] == 0)
	{
		if (MW == 0) return;
		FCFTab *F;
		if ((F = MW->getFCFTable(FCFT->item(r, 1)->text(), Mol)) == 0) return;
		addFCFTable(F, r);
	}
	if (c==1) fcf[r]->setFileName(FCFT->item(r, c)->text());
	else if (c==2) fcf[r]->setSource(FCFT->item(r, c)->text());
}

void Transition::lineTableChanged(int r, int c)
{
	if (r >= nLines || TBlock) return;
	if (Lines[r] == 0)
	{
		if (MW == 0) return;
		LineTable *L;
		if ((L = MW->getLineTable(lTable->item(r, 1)->text(), Mol)) == 0) return;
		addLineTable(L, r);
	}
	if (c==1) Lines[r]->setFileName(lTable->item(r, c)->text());
	else if (c==2) Lines[r]->setSource(lTable->item(r, c)->text());
}

void Transition::fcfTableClicked(int r, int c)
{
	if (r != nFCF || c!=0) return;
	if (FCFTB == 0) FCFTB = new ComboBox*[10];
	if (fcf == 0)
	{
		int n;
		fcf = new FCFTab*[10];
		for (n=0; n<10; n++) fcf[n] = 0;
	}
	FCFT->setCellWidget(r, c, FCFTB[nFCF] = new ComboBox(0));
	FCFTB[nFCF]->setNumber(nFCF);
	connect(FCFTB[nFCF++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkFCFTables()));
	refreshFCFBox();
}

void Transition::lineTableClicked(int r, int c)
{
	if (r != nLines || c!=0) return;
	if (LineTB == 0) LineTB = new ComboBox*[10];
	if (Lines == 0)
	{
		int n;
		Lines = new LineTable*[10];
		for (n=0; n<10; n++) Lines[n] = 0;
	}
	lTable->setCellWidget(r, c, LineTB[nLines] = new ComboBox(0));
	LineTB[nLines]->setNumber(nLines);
	connect(LineTB[nLines++], SIGNAL(currentIndexChanged(int)), this, SLOT(checkLineTables()));
	refreshLineBox();
}

void Transition::loadFCFTable(int Number)
{
	if (getFCFTable(Number) != 0) fcf[Number]->show();
}

void Transition::loadLineTable(int i)
{
	if (getLineTable(i) != 0) Lines[i]->show();
}

void Transition::lSBChanged(int i)
{
	if (lowerState != 0 ? lowerState->getStateNum() == i : i == Mol->getNumStates()) return;
	//printf("lSBChanged: StateNum=%d, i=%d\n", lowerState->getStateNum(), i);
	if (lowerState != 0) 
	{
		if (LSB != 0) disconnect(LSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		disconnect(lSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		disconnect(lowerState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	lowerState = Mol->getStateP(i);
	if (LSB != 0) LSB->setCurrentIndex(i);
	if (lowerState != 0)
	{
		if (LSB != 0) connect(LSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		connect(lSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		connect(lowerState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
    setWindowTitle("Transition " + (upperState != 0 ? upperState->getName() : QString("unknown"))
                    + " <-> " + (lowerState != 0 ? lowerState->getName() : QString("unknown")));
    Mol->updateTransitionName(transNum, i, (upperState != 0 ? upperState->getStateNum() : -1));
    updateName();
}

void Transition::LSBChanged(int i)
{
	if (lowerState != 0 ? lowerState->getStateNum() == i : i == Mol->getNumStates()) return;
	//printf("LSBChanged: StateNum=%d, i=%d\n", lowerState->getStateNum(), i);
	if (lowerState != 0) 
	{
		disconnect(LSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		disconnect(lSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		disconnect(lowerState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	lowerState = Mol->getStateP(i);
	lSB->setCurrentIndex(i);
	if (lowerState != 0)
	{
		connect(LSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		connect(lSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		connect(lowerState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
    setWindowTitle("Transition " + (upperState != 0 ? upperState->getName() : QString("unknown"))
                    + " <-> " + (lowerState != 0 ? lowerState->getName() : QString("unknown")));
    Mol->updateTransitionName(transNum, i, (upperState != 0 ? upperState->getStateNum() : -1));
    updateName();
}

void Transition::mergeLineTables()
{
	printf("Transition::mergeLineTables\n");
	QList<QTableWidgetSelectionRange> SR = lTable->selectedRanges();
	bool DelR[nLines], all = false;
	int r, l, rs, re, mT = -1, C = SR.count(), NR, NC;
	QString **data;
	for (r=0; r < nLines; r++) DelR[r] = true;
	if (C == 0)
	{
		all = true;
		mT = mLine;
	}
	else if (C == 1 && SR[0].bottomRow() == SR[0].topRow())
	{
		all = true;
		mT = SR[0].topRow();
	}
	else 
	{
		for (r=0; r < C; r++) for (l = SR[r].topRow(); l < SR[r].bottomRow(); l++)
					if (l == mLine) mT = mLine;
		if (mT == -1) mT = SR[0].topRow();
	}
	if (Lines[mT] == 0) Lines[mT] = (MW != 0 ? MW->getLineTable(lTable->item(mT, 1)->text(), Mol) : 0);
	if (Lines[mT] == 0) return;
	printf("Vor Schleife\n");
	for (r=0; r<C || all; r++) 
	{
		printf("r=%d, C=%d\n", r, C);
		if (all)
		{
			rs = 0;
			re = nLines - 1;
			all = false;
		}
		else
		{
			rs = SR[r].topRow();
			re = SR[r].bottomRow();
		}
		for (l = rs; l <= re && l < nLines; l++)
		{
			printf("l=%d, re=%d\n", l, re);
			if (Lines[l] == 0) 
				Lines[l] = (MW != 0 ? MW->getLineTable(lTable->item(l, 1)->text(), Mol) : 0);
			if (l == mT || Lines[l] == 0) continue;
			printf("Vor getData\n");
			data = Lines[l]->getData(NR = 0, NC = 0);
			printf("Vor addData\n");
			Lines[mT]->addData(data, NR, NC);
			DelR[l] = false;
			printf("Ende Schleife\n");
		}
	}
	printf("Vor Delete\n");
	deleteLineTables(DelR);
	printf("Ende mergeLineTables\n");
}

void Transition::refreshFCFBox()
{
	TBlock = true;
	int i, j, N=0;
	QString T;
	FCFTab *F;
	if (MW != 0) N = MW->getNumFCFTables();
	if (FCFB != 0)
	{
		FCFB->clear();
		for (i=0; i < nFCF; i++) FCFB->addItem(FCFTB[i]->currentText());
		FCFB->setCurrentIndex(mFCF);
	}
	for (i=0; i < nFCF; i++)
	{
		T = FCFTB[i]->currentText();
		FCFTB[i]->clear();
		FCFTB[i]->addItem(T);
		if (!T.isEmpty()) FCFTB[i]->addItem("");
	}
	for (i=0; i<N; i++) if (!(F = MW->getFCFTable(i))->isAssigned())
	{
		T = F->getName();
		for (j=0; j < nFCF; j++) FCFTB[j]->addItem(T);
		if (FCFB != 0) FCFB->addItem(T);
	}
	if (nFCF == 0 && FCFB != 0)
	{
		FCFB->addItem("");
		FCFB->setCurrentIndex(FCFB->count() - 1);
	}
	TBlock = false;
}

void Transition::refreshLineBox()
{
	TBlock = true;
	int i, j, N=0;
	QString T;
	LineTable *L;
	if (MW != 0) N = MW->getNumLineTables();
	if (LineB != 0)
	{
		LineB->clear();
		for (i=0; i < nLines; i++) LineB->addItem(LineTB[i]->currentText());
		LineB->setCurrentIndex(mLine);
	}
	for (i=0; i < nLines; i++)
	{
		T = LineTB[i]->currentText();
		LineTB[i]->clear();
		LineTB[i]->addItem(T);
		if (!T.isEmpty()) LineTB[i]->addItem("");
	}
	for (i=0; i<N; i++) if (!(L = MW->getLineTable(i))->isAssigned())
	{
		T = L->getName();
		for (j=0; j < nLines; j++) LineTB[j]->addItem(T);
		if (LineB != 0) LineB->addItem(T);
	}
	if (nLines == 0 && LineB != 0)
	{
		LineB->addItem("");
		LineB->setCurrentIndex(LineB->count() - 1);
	}
	TBlock = false;
}

void Transition::refreshStateBox()
{
    //printf("Transition::refreshStateBox()\n");
	if (Mol == 0) 
	{
		printf("Error: Mol==0!!!\n");
		return;
	}
	int n, N = Mol->getNumStates();
	int l = (lowerState != 0 ? lowerState->getStateNum() : N);
	int u = (upperState != 0 ? upperState->getStateNum() : N);
    //printf("l=%d, u=%d\n", l, u);
	QString SN;
	//printf("vor Box\n");
    lSB->blockSignals(true);
    uSB->blockSignals(true);
    if (LSB != 0) LSB->blockSignals(true);
    if (USB != 0) USB->blockSignals(true);
    lSB->clear();
	//printf("vor LSB, LSB=%d\n", LSB);
	if (LSB != 0) LSB->clear();
	//printf("vor uSB\n");
	uSB->clear();
	if (USB != 0) USB->clear();
	//printf("Vor Schleife\n");
	for (n=0; n<N; n++)
	{
		lSB->addItem(SN = Mol->getState(n));
		if (LSB != 0) LSB->addItem(SN);
		uSB->addItem(SN);
		if (USB != 0) USB->addItem(SN);
	}
	//printf("Nach Schleife\n");
	lSB->addItem("unknown");
	lSB->setCurrentIndex(l);
	if (LSB != 0) 
	{
		LSB->addItem("unknown");
		LSB->setCurrentIndex(l);
	}
	uSB->addItem("unknown");
	if (u<0 || u >= uSB->count()) 
	    printf("ACHTUNG Transition::refreshStateBox() index out of range, u=%d, uSB->count()=%d", u, uSB->count());
	else uSB->setCurrentIndex(u);
	if (USB != 0) 
	{
		USB->addItem("unknown");
		USB->setCurrentIndex(u);
	}
    lSB->blockSignals(false);
    uSB->blockSignals(false);
    if (LSB != 0) LSB->blockSignals(false);
    if (USB != 0) USB->blockSignals(false);
    //printf("Ende refreshStateBox\n");
}

void Transition::removeFCFTable(int j)
{
	bool Del[nFCF];
	int i;
	for (i=0; i < nFCF; i++) Del[i] = true;
	Del[j] = false;
	deleteFCFTables(Del);
}

void Transition::removeLineTable(int j)
{
	bool Del[nLines];
	int i;
	for (i=0; i < nLines; i++) Del[i] = true;
	Del[j] = false;
	deleteLineTables(Del);
}

void Transition::setFCFBox(ComboBox* Box)
{
	FCFB = Box;
	connect(Box, SIGNAL(activated(QString)), this, SLOT(updateFCFTable(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showFCFTable()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showLineMenu(QPoint)));
}

void Transition::setLineBox(ComboBox *Box)
{
	LineB = Box;
	connect(Box, SIGNAL(activated(QString)), this,
			 SLOT(updateLineTable(QString)));
	connect(Box, SIGNAL(DoubleClicked(int)), this, SLOT(showLineTable()));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showLineMenu(QPoint)));
}

void Transition::setLowerState(ElState *S)
{
	//printf("setLowerState S=%d\n", S);
	int n = (S != 0 ? S->getStateNum() : Mol->getNumStates());
	if (lowerState != 0) 
	{
		if (LSB != 0) disconnect(LSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		disconnect(lSB, SIGNAL(DoubleClicked(int)), lowerState, SLOT(show()));
		//disconnect(lowerState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	if (S != 0)
	{
		if (LSB != 0) connect(LSB, SIGNAL(DoubleClicked(int)), S, SLOT(show()));
		connect(lSB, SIGNAL(DoubleClicked(int)), S, SLOT(show()));
		//connect(S, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	lowerState = S;
	if (lSB->count() > n)
	{
		lSB->setCurrentIndex(n);
		if (LSB != 0) LSB->setCurrentIndex(n);
	}
    updateName();
	//printf("Ende setLowerState\n");
}

void Transition::setLSB(ComboBox *Box)
{
	LSB = Box;
	connect(Box, SIGNAL(currentIndexChanged(int)), this,
			 SLOT(LSBChanged(int)));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showLSMenu(QPoint)));
}

void Transition::setMainFCFTable(FCFTab* F)
{
	mFCF = nFCF;
	addFCFTable(F);
}

void Transition::setMainFCFTable(int BoxIndex)
{
	mFCF = BoxIndex;
	if (Mol != 0) Mol->Changed();
}

bool Transition::setMainFCFTable(QString FileName)
{
	int i;
	for (i=0; (i < nFCF ? FCFT->item(i, 1)->text() != FileName : false); i++) ;
	if (i == nFCF) return false;
	mFCF = i;
	return true;
}

void Transition::setMainLineTable(LineTable *L)
{
	mLine = nLines;
	addLineTable(L);
}

void Transition::setMainLineTable(int i)
{
	mLine = i;
	if (Mol != 0) Mol->Changed();
}

bool Transition::setMainLineTable(QString FileName)
{
	int i;
	for (i=0; (i < nLines ? lTable->item(i, 1)->text() != FileName : false); i++) ;
	if (i == nLines) return false;
	mLine = i;
	return true;
}

void Transition::setMW(MainWindow *mw, Molecule *molecule)
{
	MW = mw;
	Mol = molecule;
	ShowAction = new QAction("Show transition data", Mol);
	ShowAction->setStatusTip("Show the available data for this transition");
	connect(ShowAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowMenu = new QMenu(Mol);
	ShowMenu->addAction(ShowAction);
	ShowFCFAction = new QAction("Show FCF table", Mol);
	ShowFCFAction->setStatusTip("Show this table of Franck-Condon factors");
	connect(ShowFCFAction, SIGNAL(triggered()), this, SLOT(showFCFTable()));
	ShowFCFTAction = new QAction("Show transition data", Mol);
	ShowFCFTAction->setStatusTip("Show the available data for this transition");
	connect(ShowFCFTAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowFCFMenu = new QMenu(Mol);
	ShowFCFMenu->addAction(ShowFCFAction);
	ShowFCFMenu->addAction(ShowFCFTAction);
	ShowLineAction = new QAction("Show line table", Mol);
	ShowLineAction->setStatusTip("Show this table of transition frequencies");
	connect(ShowLineAction, SIGNAL(triggered()), this, SLOT(showLineTable()));
	ShowLineTAction = new QAction("Show transition data", Mol);
	ShowLineTAction->setStatusTip("Show the available data for this transition");
	connect(ShowLineTAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowLineMenu = new QMenu(Mol);
	ShowLineMenu->addAction(ShowLineAction);
	ShowLineMenu->addAction(ShowLineTAction);
	ShowLSSAction = new QAction("Show state properties", Mol);
	ShowLSSAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowLSSAction, SIGNAL(triggered()), this, SLOT(showLState()));
	ShowLSTAction = new QAction("Show transition data", Mol);
	ShowLSTAction->setStatusTip("Show the available data for this transition");
	connect(ShowLSTAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowLSMenu = new QMenu(Mol);
	ShowLSMenu->addAction(ShowLSSAction);
	ShowLSMenu->addAction(ShowLSTAction);
	ShowUSSAction = new QAction("Show state properties", Mol);
	ShowUSSAction->setStatusTip("Show the properties of this electronic state");
	connect(ShowUSSAction, SIGNAL(triggered()), this, SLOT(showUState()));
	ShowUSTAction = new QAction("Show transition data", Mol);
	ShowUSTAction->setStatusTip("Show the available data for this transition");
	connect(ShowUSTAction, SIGNAL(triggered()), this, SLOT(show()));
	ShowUSMenu = new QMenu(Mol);
	ShowUSMenu->addAction(ShowUSSAction);
	ShowUSMenu->addAction(ShowUSTAction);
	connect(MW, SIGNAL(LineTablesChanged()), this, SLOT(refreshLineBox()));
	connect(MW, SIGNAL(FCFTablesChanged()), this, SLOT(refreshFCFBox()));
}

void Transition::setTransNum(int N)
{
	transNum = N;
}

void Transition::setUpperState(ElState *S)
{
	//printf("setUpperState %d\n", S);
	int n = (S != 0 ? S->getStateNum() : Mol->getNumStates());
	//printf("n=%d\n", n);
	if (upperState != 0) 
	{
		if (USB != 0) disconnect(USB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		disconnect(uSB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		//disconnect(upperState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	if (S != 0)
	{
		if (USB != 0) connect(USB, SIGNAL(DoubleClicked(int)), S, SLOT(show()));
		connect(uSB, SIGNAL(DoubleClicked(int)), S, SLOT(show()));
		//connect(S, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	upperState = S;
	if (uSB->count() > n)
	{
		uSB->setCurrentIndex(n);
		if (USB != 0) USB->setCurrentIndex(n);
	}
    updateName();
	//printf("Ende setUpperState\n");
}

void Transition::setUSB(ComboBox *Box)
{
	USB = Box;
	connect(Box, SIGNAL(currentIndexChanged(int)), this,
			 SLOT(USBChanged(int)));
	connect(Box, SIGNAL(RightClicked(QPoint)), this, SLOT(showUSMenu(QPoint)));
}

void Transition::showFCFMenu(QPoint P)
{
	sFCF = mFCF;
	if (ShowFCFMenu != 0) ShowFCFMenu->popup(P);
}

void Transition::showFCFTable()
{
	if (sFCF >= 0 && sFCF < nFCF)
	{
		if (fcf[sFCF] == 0) loadFCFTable(sFCF);
		else fcf[sFCF]->show();
	}
}

void Transition::showLineMenu(QPoint P)
{
	sLine = mLine;
	if (ShowLineMenu != 0) ShowLineMenu->popup(P);
}

void Transition::showLineTable()
{
	if (sLine >= 0 && sLine < nLines)
	{
		if (Lines[sLine] == 0) loadLineTable(sLine);
		else Lines[sLine]->show();
	}
}

void Transition::showLSMenu(QPoint P)
{
	if (ShowLSMenu != 0) ShowLSMenu->popup(P);
}

void Transition::showLState()
{
	if (lowerState != 0) lowerState->show();
}

void Transition::showMenu(QPoint P)
{
	if (ShowMenu != 0) ShowMenu->popup(P);
}

void Transition::showUSMenu(QPoint P)
{
	if (ShowUSMenu != 0) ShowUSMenu->popup(P);
}

void Transition::showUState()
{
	if (upperState != 0) upperState->show();
}

void Transition::updateFCFTable(QString Name)
{
	if (MW == 0)
	{
		printf("Transition::updateFCFTable: error MW = 0!\n");
		return;
	}
	int i = FCFB->currentIndex();
	if (i < nFCF) setMainFCFTable(i);
	else
	{
		int N = MW->getNumFCFTables();
		FCFTab *F = 0;
		for (i=0; (i<N ? (F = MW->getFCFTable(i))->getName() != Name : false); i++) ;
		if (i<N) setMainFCFTable(F);
	}
	emit FCFTableChanged();
}

void Transition::updateLineTable(QString Name)
{
	if (MW == 0) 
	{
		printf("Transition::updateLineTable: error MW = 0!\n");
		return;
	}
	int i = LineB->currentIndex();
	if (i < nLines) setMainLineTable(i);
	else
	{
		int N = MW->getNumLineTables();
		LineTable *L = 0;
		for (i=0; (i<N ? (L = MW->getLineTable(i))->getName() != Name : false); i++) ;
		if (i<N) setMainLineTable(L);
	}
	emit LineTableChanged();
}

void Transition::updateName()
{
    QString Name((upperState != 0 ? upperState->getName() : "unknown") + QString(" <-> ") + (lowerState != 0 ? lowerState->getName() : "unknown"));
    setName(Name);
    setWindowTitle("Transition " + Name);
}

void Transition::uSBChanged(int i)
{
	//printf("uSBChanged\n");
	if (upperState != 0 ? upperState->getStateNum() == i : i == Mol->getNumStates()) return;
	//printf("uSBChanged: StateNum=%d, i=%d\n", upperState->getStateNum(), i);
	if (upperState != 0) 
	{
		if (USB != 0) disconnect(USB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		disconnect(uSB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		//disconnect(upperState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	upperState = Mol->getStateP(i);
	if (USB != 0) USB->setCurrentIndex(i);
	if (upperState != 0)
	{
		if (USB != 0) connect(USB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		connect(uSB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		//connect(upperState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
    setWindowTitle("Transition " + (upperState != 0 ? upperState->getName() : QString("unknown"))
                    + " <-> " + (lowerState != 0 ? lowerState->getName() : QString("unknown")));
    Mol->updateTransitionName(transNum, (lowerState != 0 ? lowerState->getStateNum() : -1), i);
    updateName();
}

void Transition::USBChanged(int i)
{
	//printf("USBChanged\n");
	if (upperState != 0 ? upperState->getStateNum() == i : i == Mol->getNumStates()) return;
	//printf("USBChanged: StateNum=%d, i=%d\n", upperState->getStateNum(), i);
	if (upperState != 0) 
	{
		disconnect(USB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		disconnect(uSB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		//disconnect(upperState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
	upperState = Mol->getStateP(i);
	if (i<0 || i >= uSB->count()) 
	    printf("ACHTUNG: Transition::USBChanged(%d) index out of range, uSB->count()=%d", i, uSB->count());
	uSB->setCurrentIndex(i);
	if (upperState != 0)
	{
		connect(USB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		connect(uSB, SIGNAL(DoubleClicked(int)), upperState, SLOT(show()));
		//connect(upperState, SIGNAL(PotentialChanged()), this, SLOT(updateFCF()));
	}
    setWindowTitle("Transition " + (upperState != 0 ? upperState->getName() : QString("unknown"))
                    + " <-> " + (lowerState != 0 ? lowerState->getName() : QString("unknown")));
    Mol->updateTransitionName(transNum, (lowerState != 0 ? lowerState->getStateNum() : -1), i);
    updateName();
}

bool Transition::writeData(QString FileName)
{
    if (Mol == 0) return MDIChild::writeData(FileName);
	if (QMessageBox::information(this, "MolSpektAnalysis", "The transition is part of the molecule " + Mol->getName() 
			+ " and connot be saved separately. Instead the data of the molecule will be written to disk.", 
			QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok) == QMessageBox::Ok)
		return Mol->writeData(FileName);
	return true;
}
