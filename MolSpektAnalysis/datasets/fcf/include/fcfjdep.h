//
// C++ Interface: FCFJDep
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2010 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FCFJDEP_H
#define FCFJDEP_H


#include <QWidget>

class DiagWindow;
class Potential;
class Molecule;
class MainWindow;

class QPushButton;
class QComboBox;
class QCheckBox;


class FCFJDep : public QWidget
{
	Q_OBJECT
	
	public:
    FCFJDep(MainWindow *MW);
	
	private slots:
		void addToPlot();
		void newPlot();
		void clearPlot();
		void molChanged(int i);
		void lPotChanged();
		void uPotChanged();
		void vChanged();
		void minJChanged(int i);
		void maxJChanged(int i);
		
	private:
		DiagWindow *lDiag;
		QPushButton *addB, *newB, *clearB, *closeB;
		QComboBox *MolB, *IsoB, *lPotB, *uPotB, *JssB, *vsB, *vssB, *minJB, *maxJB;
		QCheckBox *Fit;
		Potential *lPot, *uPot;
		double ***lData, ***uData;
		int Nvss, Nvs, NJss, NJs, NIso, aMJ, amJ;
		Molecule *Mol;
		MainWindow *mw;
};

#endif // FCFJDEP_H
