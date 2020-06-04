//
// C++ Interface: Lorentzian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LORENTZIAN_H
#define LORENTZIAN_H


#include "fitobject.h"


class Lorentzian : public FitObject
{
public:
    Lorentzian();
	Lorentzian(double *x, double *y, double *Sig, int N);
	double getE();
    void setData(double* x, double* y, double* Sig, int N);
	
protected:
	void getDerivatives(double **deriv);
    bool getCalcY(double *YCalc);
	void getPar(double *Par);
	void setPar(double *Par);
	void updatePar(double *C);
	
private:
	double A, E, Gamma, AWF, EWF, GWF;
};

#endif
