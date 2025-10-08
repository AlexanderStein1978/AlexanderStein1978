//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FELDIALOG_H
#define FELDIALOG_H


#include <QDialog>

class MainWindow;
class ElState;
class Molecule;

class QComboBox;
class QLineEdit;



class FELDialog : public QDialog
{
	Q_OBJECT
	
public:
	FELDialog(MainWindow *MW, double initMinH);
	Molecule *getMol();
	ElState *getState();
	int getIso();
	int getMLPB();
	double getTol();
	double getMinH();
	
	static bool DataA(MainWindow *MW);
	static bool DataA(Molecule *Mol);
	static bool DataA(Molecule *Mol, ElState *State);
	static bool DataA(Molecule *Mol, ElState *State, int Iso);
	
private slots:
	void MolChanged(QString Name);
	void StateChanged(QString Name);
	
private:
	Molecule *Mol;
	ElState *State;
	MainWindow *MW;
	QComboBox *MolBox, *StateBox, *IsoBox;
	QLineEdit *MLPBEdit, *TolEdit, *MinHEdit;
	QPushButton *OK, *Cancel;
};

#endif
