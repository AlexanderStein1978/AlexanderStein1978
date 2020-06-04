//
// C++ Interface: atom
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ATOM_H
#define ATOM_H

#include "constants.h"

#include "mdichild.h"

class QString;
class QLineEdit;
class QTableWidget;
class QPushButton;

class MainWindow;

struct Isotop 
{
	int nNuc;
	QString Mass;
	QString NA;
	QString CoreSpin;
};

class Atom : public MDIChild
{
	Q_OBJECT

public:
	Atom(MainWindow *MW = 0);
    ~Atom();
	bool readData(QString Filename);
	QString getChSymb();
	QString getSourceData();
	int getNP();
	int getnIso();
	int getnNuc(int iIso);
	double getIsoMass(int iIso);
	double getIsoNA(int iIso);
	float getCoreSpin(int iIso);
public slots:
	bool writeData(QString Filename = "");
private slots:
	void Reset();
	void OK();
signals:
    void closeThis();
private:
	QString ChSymb, SourceData;
	int NP, nIso;
	Isotop Iso[MaxIso];
	QLineEdit *SourceIn, *NameIn, *SymIn, *npIn;
	QTableWidget *IsoTab;
	QPushButton *SaveB, *OKB, *ResetB, *CloseB;
};

#endif
