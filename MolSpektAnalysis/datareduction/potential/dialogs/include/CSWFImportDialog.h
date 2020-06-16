//
// C++ Interface: CSWFImportDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CSWFIMPORTDIALOG_H
#define CSWFIMPORTDIALOG_H


#include <QDialog>

class MainWindow;
class Molecule;
class ElState;
class Molecule;

class QListWidget;
class QComboBox;
class QLineEdit;


class CSWFImportDialog : public QDialog
{
	Q_OBJECT
	
public:
    CSWFImportDialog(MainWindow *MW, QString InitDir);
	void getData(Molecule *&Mol, ElState **&States, int *&Components, int &NStates, int *&Iso, 
				 QString *&IsoDirs, int &NIso, QString &PotentialFile, QString &T5_E_BvFile);

private slots:
	void molChanged(int i);
	void addState();
	void addIso();
	void removeState();
	void removeIso();
	void selectPotFile();
	void selectT5_E_Bv();
	void selectIsoDir();
	
private:
	
	int *ST, *CT;
	MainWindow *mw;
	Molecule *mol;
	QListWidget *SList, *IsoList;
	QComboBox *SBox, *MolBox, *IsoBox;
	QPushButton *AddState, *RemoveState, *AddIso, *RemoveIso, *SelectPotFile, *SelectT5_E_Bv, *SelectIsoDir, *OK, *Cancel;
	QLineEdit *PotFile, *T5_E_Bv, *IsoDir;
	QList<int> IL1, IL2, SL1, SL2;
};

#endif
