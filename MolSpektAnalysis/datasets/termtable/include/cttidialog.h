//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CTTIDIALOG_H
#define CTTIDIALOG_H


#include <QDialog>


class QListWidget;
class QComboBox;
class QPushButton;
class QLineEdit;
class QRadioButton;

class MainWindow;
class Molecule;
class ElState;


class CTTIDialog : public QDialog
{
	Q_OBJECT
	
public:
	CTTIDialog(MainWindow *MW);
	void getData(Molecule *&Mol, ElState **&States, int &NStates, int &Mv, int &MJ, bool &singleTable);

private slots:
	void molChanged(int i);
	void add();
	void remove();
	
private:
	
	MainWindow *mw;
	Molecule *mol;
	QListWidget *SList;
	QComboBox *SBox, *MolBox;
	QPushButton *Add, *Remove, *OK, *Cancel;
	QLineEdit *MaxJ, *Maxv;
	QRadioButton *STable, *MTables;
};

#endif
