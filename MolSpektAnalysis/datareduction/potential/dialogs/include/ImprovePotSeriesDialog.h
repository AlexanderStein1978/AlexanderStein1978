//
// C++ Interface: ImprovePotSeriesDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef IMPROVEPOTSERIESDIALOG_H
#define IMPROVEPOTSERIESDIALOG_H


#include <QDialog>

class Potential;
class FitData;
class MainWindow;
class ElState;
class Molecule;

class QComboBox;
class QLineEdit;


class ImprovePotSeriesDialog : public QDialog
{
	Q_OBJECT
	
public:
	ImprovePotSeriesDialog(MainWindow *MW, QString Dir, Potential *OPot);
	bool exec(ElState *&St, FitData *&FDat, QString &PotDir, double &Threshold, int &NumParFits);
	
private slots:
	void MolChanged(int i);
	void StateChanged(int i);
	void sPotDir();
	
private:
	MainWindow *mw;
	Molecule *Mol;
	ElState *St;
	QComboBox *MoleculeB, *StateB, *FitDataB;
	QLineEdit *ThresholdE, *NumParFitsE, *PotDirE;
};

#endif
