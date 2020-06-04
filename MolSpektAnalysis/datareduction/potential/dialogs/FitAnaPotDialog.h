//
// C++ Interface: FitAnaPotDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FITANAPOTDIALOG_H
#define FITANAPOTDIALOG_H


#include <QDialog>


class FitAnaPotDialog : public QDialog
{
	Q_OBJECT
public:
    FitAnaPotDialog(bool improve, double Ri, double Ra, double iExp, int NC, int NLRC, int *pLRC, double hCi_cP, 
					PotentialType CurrentPotential, QWidget* parent = 0);
	bool exec(bool &improve, bool &useLeveMrq, bool &useSVD, double &Ri, double &Ra, double &iExp, int &NC, int &NLRC, int *&pLRC, 
			  double &Rmin, double &Rmax, int &NumP, double &bmin, double &bmax, double &bstep, bool &adjustTEWeightFact, 
		      double &hCi_cP, PotentialType &NewType);
	
private slots:
	void improveChanged(bool State);
	void addCoeff();
	void removeCoeff();
	void typeChanged(int Index);
	
private:
	QListWidget *CoeffList;
	QComboBox *CoeffBox, *TypeBox;
	QPushButton *Add, *Remove;
	QLineEdit *RiE, *RaE, *iExpE, *NCE, *RminE, *RmaxE, *NumPE, *bminE, *bmaxE, *bstepE, *hCi_cPE;
	QCheckBox *Improve, *UseLeveMrq, *UseSVD, *AdjustTEWF;
	PotentialType CurPotential;
};

#endif
