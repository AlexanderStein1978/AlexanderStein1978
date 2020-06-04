//
// C++ Implementation: atom
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//
#include <QLabel>
#include <QLineEdit>
#include <QTableWidget>
#include <QPushButton>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QGridLayout>
#include "atom.h"

Atom::Atom(MainWindow *mw) : MDIChild(AtomData, mw, "Atoms (*.atom)", ".atom")
{
	printf("Atom::Atom\n");
	NP = nIso = 0;
	
	setMinimumSize(400, 400);

    QGridLayout *L = new QGridLayout(this);
	QLabel *Label1 = new QLabel("Source:", this);
    L->addWidget(Label1, 0, 0);
	
	SourceIn = new QLineEdit("", this);
    L->addWidget(SourceIn, 0, 1, 1, 7);
	connect(SourceIn, SIGNAL(textEdited(QString)), this, SLOT(Changed()));
	 
	QLabel *Label2 = new QLabel("Name:", this);
    L->addWidget(Label2, 1, 0);
	
	NameIn = new QLineEdit("", this);
    L->addWidget(NameIn, 1, 1, 1, 3);
	connect(NameIn, SIGNAL(textEdited(QString)), this, SLOT(Changed()));
	
	QLabel *Label3 = new QLabel("Symbol:", this);
    L->addWidget(Label3, 1, 4);
	
	SymIn = new QLineEdit("", this);
    L->addWidget(SymIn, 1, 5);
	connect(SymIn, SIGNAL(textEdited(QString)), this, SLOT(Changed()));
	
	QLabel *Label4 = new QLabel("n.p.:", this);
    L->addWidget(Label4, 1, 6);
	
	npIn = new QLineEdit("", this);
    L->addWidget(npIn, 1, 7);
	connect(npIn, SIGNAL(textEdited(QString)), this, SLOT(Changed()));
	
    L->setRowMinimumHeight(2, 20);

	QLabel *Label5 = new QLabel("Isotopes:", this);
    L->addWidget(Label5, 3, 0, 1, 2);
	
	IsoTab = new QTableWidget(MaxIso, 4, this);
    L->addWidget(IsoTab, 4, 0, 1, 8);
	IsoTab->setHorizontalHeaderLabels(QStringList() << "mass number" << "mass [u]" 
			<< "rel. na" << "core spin");
	//IsoTab->setColumnWidth(1, IsoTab->width() - IsoTab->columnWidth(0) 
	//		- IsoTab->columnWidth(2) - 43);
	connect(IsoTab, SIGNAL(cellChanged(int, int)), this, SLOT(Changed()));
	
    L->setRowMinimumHeight(5, 20);

	SaveB = new QPushButton("Save", this);
    L->addWidget(SaveB, 6, 0, 1, 2);
	connect(SaveB, SIGNAL(clicked()), this, SLOT(writeData()));
	
	OKB = new QPushButton("OK", this);
    L->addWidget(OKB, 6, 2, 1, 2);
	connect(OKB, SIGNAL(clicked()), this, SLOT(OK()));
	
	ResetB = new QPushButton("Reset", this);
    L->addWidget(ResetB, 6, 4, 1, 2);
	connect(ResetB, SIGNAL(clicked()), this, SLOT(Reset()));
	
	CloseB = new QPushButton("Close", this);
    L->addWidget(CloseB, 6, 6, 1, 2);
    connect(CloseB, SIGNAL(clicked()), this, SIGNAL(closeThis()));
	
    L->setColumnStretch(0, 1);
    L->setColumnStretch(1, 1);
    L->setColumnStretch(2, 1);
    L->setColumnStretch(3, 1);
    L->setColumnStretch(4, 1);
    L->setColumnStretch(5, 1);
    L->setColumnStretch(6, 1);
    L->setColumnStretch(7, 1);

	setWindowTitle("New Atom");
}


Atom::~Atom()
{
}

bool Atom::writeData(QString NFilename)
{
	int i;
	OK();
	QFile Datei(NFilename);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
    S << "Source of masses: " << SourceData << "\n";
	S << "Name: " << getName() << "\n";
	S << "Chemical symbol: " << ChSymb << "\n";
    S << "Atomic number: " << NP << "\n\n";
	S << "Table of isotopes:\n";
    S << "Number of nuclides | mass [u] | relative abundance | nuclear spin\n";
	for (i=0; i<nIso; i++) 
		S << Iso[i].nNuc << " | " << Iso[i].Mass << " | " << Iso[i].NA << " | " << Iso[i].CoreSpin << "\n";
	return true;
}

bool Atom::readData(QString FileName)
{
	printf("Atom::readData\n");
	int n, m;
	bool nT = true;
	nIso = 0;
	QString Buffer;
	QFile Datei(FileName);
	if (!read(&Datei)) return false;
	QTextStream S(&Datei);
	while (!S.atEnd())
	{
		Buffer = S.readLine();
		if (nT)
		{
			if (Buffer.indexOf("Source", 0, Qt::CaseInsensitive) != -1)
			{
				for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
				SourceData = Buffer.right(Buffer.length() - n);
			}
			else if (Buffer.indexOf("Name", 0, Qt::CaseInsensitive) != -1)
			{
				for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
				setName(Buffer.right(Buffer.length() - n));
			}
			else if (Buffer.indexOf("symbol", 0, Qt::CaseInsensitive) != -1)
			{
				for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
				ChSymb = Buffer.right(Buffer.length() - n);
			}
            else if (Buffer.indexOf("Atomic number:", 0, Qt::CaseInsensitive) != -1 || Buffer.indexOf("proton", 0, Qt::CaseInsensitive) != -1)
				NP = Buffer.right(Buffer.length() - Buffer.indexOf(":") - 1).toInt();
			else if (Buffer.indexOf("Table", 0, Qt::CaseInsensitive) != -1) nT = false;
		}
		else 
		{
			n = Buffer.indexOf("|");
			if ((m = Buffer.left(n).toInt()) > 0 && nIso < MaxIso)
			{
				Iso[nIso].nNuc = m;
				for (m = Buffer.indexOf("|", n + 1); Buffer[m].isSpace() || Buffer[m] == '|'; m--) ;
				while (Buffer[n].isSpace() || Buffer[n] == '|') n++;
				Iso[nIso].Mass = Buffer.mid(n, m - n + 1);
				for (m++; Buffer[m].isSpace() || Buffer[m] == '|'; m++) ;
				for (n = Buffer.indexOf("|", m + 1); 
								 (n != -1 ? Buffer[n].isSpace() || Buffer[n] == '|' : false); n--) ;
				if (n == -1) n = Buffer.length();
				Iso[nIso].NA = Buffer.mid(m, n - m + 1);
				if (n < (m = Buffer.length())) 
				{
					for (n++; Buffer[n].isSpace() || Buffer[n] == '|'; n++) ;
					Iso[nIso].CoreSpin = Buffer.right(m - n);
				}
				nIso++;
				//printf("Buffer = %s\n", Buffer.ascii());
			}
		}
	}
	Reset();
	Saved();
	return true;
}

void Atom::OK()
{
	QTableWidgetItem *F;
	if (getName() != NameIn->text()) setName(NameIn->text());
	ChSymb = SymIn->text();
	SourceData = SourceIn->text();
	NP = npIn->text().toInt();
	for (nIso = 0; (nIso < MaxIso ? (0 != (F = IsoTab->item(nIso, 0)) ? !F->text().isEmpty() : false) : false);
			       nIso++)
	{
		Iso[nIso].nNuc = F->text().toInt();
		Iso[nIso].Mass = (0 != (F = IsoTab->item(nIso, 1)) ? F->text() : "0.0");
		Iso[nIso].NA = (0 != (F = IsoTab->item(nIso, 2)) ? F->text() : "0.0");
		Iso[nIso].CoreSpin = (0 != (F = IsoTab->item(nIso, 3)) ? F->text() : "-1");
	}
}

void Atom::Reset()
{
	int i, j;
	QTableWidgetItem *F;
	NameIn->setText(getName());
	SymIn->setText(ChSymb);
	SourceIn->setText(SourceData);
	npIn->setText(QString::number(NP));
	for (i = 0; i < nIso; i++)
	{
		IsoTab->setItem(i, 0, new QTableWidgetItem(QString::number(Iso[i].nNuc)));
		IsoTab->setItem(i, 1, new QTableWidgetItem(Iso[i].Mass));
		IsoTab->setItem(i, 2, new QTableWidgetItem(Iso[i].NA));
		IsoTab->setItem(i, 3, new QTableWidgetItem(Iso[i].CoreSpin));
	}
	for (i = nIso; i < MaxIso; i++) for (j = 0; j < 3; j++) 
			if ((F = IsoTab->item(i, j)) != 0) F->setText("");
}

QString Atom::getChSymb()
{
	return ChSymb;
}

float Atom::getCoreSpin(int iIso)
{
	if (iIso < 0 && iIso >= nIso)
	{
		printf("Atom::getnNuc: error! iIso out of valid range!\n");
		return 0;
	}
	int n;
	if ((n = Iso[iIso].CoreSpin.indexOf("/2")) == -1) return Iso[iIso].CoreSpin.toFloat();
	return Iso[iIso].CoreSpin.left(n).toFloat() / 2;
}

QString Atom::getSourceData()
{
	return SourceData;
}

int Atom::getNP()
{
	return NP;
}

int Atom::getnIso()
{
	return nIso;
}

int Atom::getnNuc(int iIso)
{
	if (iIso >= 0 && iIso < nIso) return Iso[iIso].nNuc;
	printf("Atom::getnNuc: error! iIso out of valid range!\n");
	return 0;
}

double Atom::getIsoMass(int iIso)
{
	if (iIso >= 0 && iIso < nIso) return Iso[iIso].Mass.toDouble();
	printf("Atom::getIsoMass: error! iIso out of valid range!\n");
	return 0.0;
}

double Atom::getIsoNA(int iIso)
{
	if (iIso >=0 && iIso < nIso) return Iso[iIso].NA.toDouble();
	printf("Atom::getIsoNA: error! iIso out of valid range!\n");
	return 0.0;
}
