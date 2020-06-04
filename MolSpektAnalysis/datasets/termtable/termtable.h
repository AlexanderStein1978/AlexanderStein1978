//
// C++ Interface: termtable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TERMTABLE_H
#define TERMTABLE_H

#include "tablewindow.h"

#include <QDialog>


class Molecule;
class IsoTab;

class QListWidget;
class QRadioButton;

struct Perturbation;


class TermTable : public TableWindow
{
public:
    TermTable(MainWindow *MW = 0);
    ~TermTable();
	bool readData(QString Filename);
	bool writeData(QString Filename = "");
	double ****getData();
	int getNumComp();
	int getNumIso();
	int getMaxv();
	int getMaxJ();
	ElState *getElState();
	void setMolecule(Molecule *mol, ElState *nState);
	void setData(double ****nData, int numComp, int numIso, int maxv, int maxJ, int *CompZ = 0, int nStates = 0, double *****MixC = 0);
	void GetIsoZ(int IsoI, int &NIso1, int &NIso2);
	int *getIsoT();
	void getCompT(int &Max, int *&Trans);
	void getThermPopulation(double **&Pop, int &nv, int &nJ, double T, int Iso, int Comp);
	int getNPerturbations(int c = -1, int Iso = -1, int v = -1);
	Perturbation *getPerturbation(int c, int Iso, int v, int Number);
	Perturbation *getPerturbation(int N);
	void setPerturbations(int NPerturbations, Perturbation* NewPerturbations);
	virtual void setName(QString name);
	void getViewnE(int*& Js, double*& E, int& N);
	void getSelE(int*& Js, double*& E, int &N);
	void setViewnLevels(MDIChild *Viewer, int C, int I, int v, int *J, int N);
	void getMaxAndAverageDeviation(TermTable *Other, int comp, int iso, int vMax, 
								   int JMax, double &MaxDev, double &AvDev);
	
	inline int *getCompZ()
	{
		return tDat->getCompZ();
	}
	
	inline int *getIsoZ()
	{
		return tDat->getIsoZ();
	}
	
	inline void setIsoZ(int *Z)
	{
		tDat->setIsoZ(Z);
	}
	
	inline int getNumLines()
	{
		return tDat->rowCount();
	}
	
	inline int getNumStates()
	{
		return tDat->getnumStates();
	}
	
	inline double*****getMixCoefficients()
	{
		return tDat->getMixCoeff();
	}
	
protected:
	TermData *tDat;
	
private:
	bool OpenXnfitTO(QString Dateiname);
	void AssignIso();
	ElState *State;
	int NPerturbations;
	Perturbation *Perturbations;
};

#endif
