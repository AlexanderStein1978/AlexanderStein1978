//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef DUNTABLE_H
#define DUNTABLE_H

#include "constants.h"

#include "tablewindow.h"
#include "elstate.h"

class IsoTab;
class TableLine;
class TermTable;
class Transition;


/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/
class DunTable : public TableWindow
{
	Q_OBJECT
public:
	DunTable(MainWindow *MW = 0, Molecule *Mol = 0);
    ~DunTable();
	bool readData(QString Filename);
	bool writeData(QString Filename = "");
	void exportTF();
	void getData(int &KMax, int *&LMax, double **&Par, double **&Korr, double **&SpinRot, double **&adCorr);
	void setData(int numCoefficients, int *type, int *k, int *l, double *Coefficient);
	void setElState(ElState *ES);
	ElState *getElState();
	void Fit(int mk = 0, int ml = 0, bool RemDataP = false, int vp1 = 0, int Jp1 = 0, int vp2 = 0, 
			 int Jp2 = 0);
	void getLinePoints(int &vp1, int &Jp1, int &vp2, int &Jp2);
	void setLinePoints(int vp1, int Jp1, int vp2, int Jp2);
	int getMaxk();
	int getMaxl();
	QString getTexTable();
	int *getmaxv();
	void addTableLine();
	double getwe();
	void removeFinestructure(TableLine *&TL, int &N);
	bool isSpinRAv();
	void removeLambdaDoubling(TableLine *&TL, int &N);
	void applyAdiabaticCorrection(TableLine *&TL, int &N);
	bool isLambdaDoublingAv();
	void calcRe(double &Re, double &Err);
	void calcTeY00(double &Te, double &TeErr, double &Y00, double &Y00Err);
	void testKratzer(double &Y02f, double &Y02c);
	void getBandConstants(double **&C, double **&err, int &nC, int &NLD, int &NSR, int &NAD, int &nv);
    int getBe(double &o_Be) const;
	
	inline void AddRow()
	{
		addTableLine();
	}
	
public slots:
	void calcTermEnergies(TermTable *&TT, bool show = false, int vMax = -1, int JMax = -1);
	void itemChanged(int row, int column);
private:
	int **getDunList(int &nC, int &NC, bool &FS, int &mk, int &ml);
	void getStates(int &N, int &nP, int &eC, Transition **ES, int &gC, Transition **GS, int &mv, 
				   int &mJ);
	void calcXY(Transition **GS, int gC, Transition **ES, int eC, int NP, int xx[][4], double *yy,
				double *ssig, double *SS);
	void calcXY(FitData *FD, int NP, int &N, int xx[][5], double *yy, double *ssig, 
				int *RowN, double *SS);
	double fitDunSet(int N, int nP, FitData *FD, int **DunList, 
				    int Ml, int Mk, int &nC, double ***vF, double ***JF, double *Res, double *err, 
	   				bool RemD, IsoTab *Iso);
	double fitDunCoeff(int N, int NP, FitData *FD, int **dunList,
					   int mL, int mk, int NC, double ***vf, double ***Jf, double *Res, double *err, IsoTab *Iso);
	void eIsoO(int *IsoO, int &NObsIso, int NIso, bool *ObsIso);
	double ***cvF(int NIso, IsoTab *IsoT, int maxv, int maxk);
	double ***cJF(int NIso, IsoTab *IsoT, int maxJ, int maxl, int Omega);
	double ****cSpinRotJF(IsoTab *IsoT, int maxJ, int maxl);
//	void pLA(TableLine *&Lines, int &NLines);
	
	int numCoefficients, vp1, Jp1, vp2, Jp2;
	ElState *State;
	bool BlockIC;
};

#endif
