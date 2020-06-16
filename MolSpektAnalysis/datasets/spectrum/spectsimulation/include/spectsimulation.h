//
// C++ Interface: spectsimulation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SPECTSIMULATION_H
#define SPECTSIMULATION_H


#include "DiagWindow.h"

class Molecule;
class ElState;
class IsoTab;

struct SimLine;


class SpectSimulation : public DiagWindow
{
	Q_OBJECT
public:
    SpectSimulation(MainWindow *MW);
    ~SpectSimulation();
	double voigtP(double R, double f, double I, double RRes, double IRes, double DF, 
				  double L_gsq) const;
	void voigt(double f, double I, double &Rs, double &Re, double *&Data, double RRes, double IRes,
			   double DF, double L_gsq) const;
	void simulateAbsorption(Molecule *Mol, double Rs, double Re, double RRes, double IRes, double T,
							double *Data, int &NLines, SimLine *&Lines) const;
	void cExcProfile(double Cf, double I, double LF, double &Rs, double &Re, double *&Data, 
					 double RRes, double IRes, double DF, double L_gsq) const;
	void flLine(double Cf, double Lf, double I, double &Rs, double &Re, double *&Data, double RRes,
				double IRes, double L_gsq, double EPRs, double EPRe, double *EPData) const;
	void simulateLIF(Molecule *Mol, double Rs, double Re, double RRes, double IRes, double T, 
					 double LF, double *Data, int &NLines, SimLine *&Lines) const;
	void simulateFilteredAbsorption(Molecule *Mol, double Rs, double Re, double RRes, double IRes,
									double T, double Ws, double We, double *Data, int &NLines,
		 							SimLine *&Lines) const;
public slots:
	void RefreshMolecules();
	bool readData(QString FileName);
	bool writeData(QString FileName);
	
private slots:
	void Simulate();
	void MolBoxChanged(QString NMol);
	void SimBoxChanged(int Sim);
	
private:
	void PSpektrum(QPainter &P, const QRect &R, bool PrintFN);
	void getStates(Molecule *Mol, int &N, ElState **&ES) const;
	void heapSortLines(SimLine *&Lines, int NLines) const;
	
	Molecule *Mol;
	int NumLines, Sim;
	SimLine *Lines;
	IsoTab *Iso;
};

#endif
