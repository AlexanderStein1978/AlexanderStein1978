//
// C++ Interface: SFQSCalcDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SFQSCALCDIALOG_H
#define SFQSCALCDIALOG_H


#include <QDialog>

class MainWindow;
class ElState;
class Potential;

class QComboBox;
class QLineEdit;


class SFQSCalcDialog : public QDialog
{
	Q_OBJECT
	
public:
	SFQSCalcDialog(MainWindow *MW, QString Dir);
	void exec(ElState *&St, Potential *&OPot, QString& PotDir, QString &FDatDir, double &SFQSRad, int &NumParFits);
	
private slots:
	void MolChanged(int n);
	void StateChanged(int n);
	void sPotDir();
	void sFDatDir();
	
private:
	QComboBox *MolBox, *StateBox, *PotBox;
	QLineEdit *PotDirE, *FDatDirE, *SFQSRadE, *NumParFitsE;
	QPushButton *sPotDirB, *sFDatDirB, *OK, *Cancel;
	QString Dir;
	MainWindow *MW;
};

#endif
