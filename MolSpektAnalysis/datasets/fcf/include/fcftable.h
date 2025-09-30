//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FCFTABLE_H
#define FCFTABLE_H


#include "tablewindow.h"
#include "MainWindow.h"


class FCFTable : public TableWindow
{
	Q_OBJECT
public:
    FCFTable(MainWindow *mw, Molecule *Mol = 0, QString LState = "", QString UState = "");
    ~FCFTable();
private slots:
	void FillTable();
	void UpdateMolData();
	void UpdateJ();
	void ViewChanged(int i);
private:
	int JDiff;
};

#endif
