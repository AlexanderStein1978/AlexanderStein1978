//
// C++ Interface: CoupledSineWaveFunc
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef COUPLEDSINEWAVEFUNC_H
#define COUPLEDSINEWAVEFUNC_H


#include "mdichild.h"

class Molecule;
class ElState;
class NaturalSpline;


class CoupledSineWaveFunc : public MDIChild
{
public:
    CoupledSineWaveFunc(MainWindow *MW, Molecule *Mol);
    ~CoupledSineWaveFunc();
	bool readData(QString FN = "");
	bool writeData(QString FN = "");
	int getNChannels();
	int getNLevels(int I, int J);
	void getStateComp(int Channel, ElState *&State, int &Comp);
	bool getWaveFunc(int I, int J, int v, int Channel, double RStart, double RStop, int NPoints, double *Data, double &E, 
					 int &MMC);
    bool getWaveFuncs(int I, int J, int Channel, double RStart, double RStop, int NPoints, double **WF, double *E);
    bool getWaveFuncs(const int I, const int J, int &Nv, const int NChannels, const double RStart, const double RStop, const int NPoints,
                      double ***&WFs, int*& ri, int*& ra, int*& P, double *&F, double *&N);
	bool getWaveFuncs(int I, int v, int *Channels, int NChannels, int NJ, int *J, double RStart, 
					  double RStop, int NPoints, double ***WF, double *E, int *vAssignment,
					  int *ChannelAssignment);
	bool getWaveFuncs(int I, int J, double RStart, double RStop, int NPoints, double ***WF, double *E);
	void getv(int I, int J, int &N, int *&v);
	int getv(int I, int J, double E);
	bool isJA(int I, int J);
	void getIso(int &N, int *&I);
	void setData(int NI, int NJ, int Nv, int NChan, int NCoeff, ElState **St, int *C, int *I, double *R, double *****Data,
				 double ****MC, double ***E);
	double ****getEnergy();
	
	inline int getNJ()
	{
		return nJ;
	}
	
	inline int getNv()
	{
		return nv;
	}
	
	inline QString getSource()
	{
		return Source;
	}
	
	inline void setSource(QString S)
	{
		Source = S;
	}
	
	inline void Assign()
	{
		ASC++;
	}
	
	inline bool Delete()
	{
		if (--ASC == 0) return true;
		return false;
	}
	
	inline int getIsoT(int I)
	{
		int n;
		for (n=0; n < nI; n++) if (I == Iso[n]) return n;
		return -1;
	}
	
private:
	void DestroyData();
	void CalcBasis(double rMin, double rMax, int NPoints);
	void CalcvCA();
	
	int nI, nJ, nv, nChan, nCoeff, *Iso, *SC, ASC, nPoints, ***vA, ***CA;
	double *R, *****Data, ****MC, ***E, **Basis, RMin, RMax;
	ElState **States;
	Molecule *Mol;
	QString Source;
	NaturalSpline *MapFuncs;
};

#endif
