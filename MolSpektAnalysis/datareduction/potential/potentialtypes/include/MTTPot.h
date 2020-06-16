//
// C++ Interface: MTTPot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MTTPOT_H
#define MTTPOT_H


#include "PotWorker.h"


class MTTPot : public PotWorker
{
public:
	MTTPot(PotFit *Fit);
	MTTPot(const MTTPot &C);
	~MTTPot();
	void getMinimum(double &R, double &T, int c=0);
	double GetDe(int c=0);
	double GetUinf(int c=0);
	MTTPot *scalePotential(double newRe, double newDe);
	void setCoefficients(double B, double alpha1, int NA, double *A, double alpha, double beta, double gamma, 
						 int NLRC, int *pLRC, double *LRC, double Uinf, int NC = 0,
						 double *C = 0, double Calpha = 0.0, double Cbeta = 0.0,
						 int NCLRC = 0, int *pCLRC = 0, double *CLRC = 0, 
					     double CInf = 0.0);
	void getCoefficients(double &B, double &alpha1, int &NA, double *&A, double &alpha, double &beta, double &gamma,
						 int &NLRC, int *&pLRC, double *&LRC, double &Uinf, int &NC,
						 double *&C, double &Calpha, double &Cbeta, int &NCLRC,
					     int *&pCLRC, double *&CLRC, double &CInf);
	PotentialData *getPotentialData();
	
	inline bool isFSA()
	{
		if (NC > 0 && Calpha != 0.0 && Cbeta != 0.0 && NCLRC > 0) return true;
		return false;
	}
	
protected:
	void getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients);
	double Point(double r, int c=0);
	
private:
	int NA, NLRC, *pLRC, NC, NCLRC, *pCLRC;
	double *A, alpha, beta, gamma, *LRC, Uinf, *C, Calpha, Cbeta, *CLRC, CInf, B, alpha1;
};

#endif
