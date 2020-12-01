//
// C++ Implementation: Lorentzian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "lorentzian.h"
#include "utils.h"
#include "tools.h"

#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>

#include <QFile>
#include <QTextStream>


using std::isnan;

Lorentzian::Lorentzian() : FitObject(3)
{
	AWF = EWF = GWF = A = E = Gamma = 0.0;
}

Lorentzian::Lorentzian(double* x, double* y, double* Sig, int N) : FitObject(3)
{
	setData(x, y, Sig, N);
}

bool Lorentzian::getCalcY(double* YCalc) const
{
    int n;
	double AG = A * Gamma;
	for (n=0; n < nData; n++) YCalc[n] = AG / (sqr(E - X[n]) + Gamma);
    return true;
}

void Lorentzian::getDerivatives(double** deriv)
{
    double ME = 0.0, MA = 0.0, MG = 0.0, D, DSq, Div, AG = A * Gamma;
	int n;
	for (n=0; n < nData; n++)
	{
		D = E - X[n];
		DSq = sqr(D);
		Div = 1.0 / (DSq + Gamma);
		deriv[n][0] = Gamma * Div;
		deriv[n][1] = -2 * AG * D * sqr(Div);
		deriv[n][2] = A * Div * (1.0 - Gamma * Div);
		if (fabs(deriv[n][0]) > MA) MA = fabs(deriv[n][0]);
		if (fabs(deriv[n][1]) > ME) ME = fabs(deriv[n][1]);
		if (fabs(deriv[n][2]) > MG) MG = fabs(deriv[n][2]);
	}
	AWF = 1.0 / MA;
	EWF = 1.0 / ME;
	GWF = 1.0 / MG;
	for (n=0; n < nData; n++)
	{
		deriv[n][0] *= AWF;
		deriv[n][1] *= EWF;
		deriv[n][2] *= GWF;
	}
}

double Lorentzian::getE()
{
	return E;
}

void Lorentzian::getPar(double* Par)
{
    Par[0] = A;
	Par[1] = E;
	Par[2] = Gamma;
}

void Lorentzian::setData(double* x, double* y, double* Sig, int N)
{
    int n, m=0;
	AWF = EWF = GWF = 0.0;
	for (n=0, A=0.0; n<N; n++) if (y[n] > A)
	{
		m=n;
		A = y[n];
	}
	E = x[m];
	if (N>2 && m>0 && m<N-1) Gamma = sqr(x[m+1] - x[m-1]);
	else Gamma = 1;
	FitObject::setData(x, y, Sig, N);
}

void Lorentzian::setPar(double* Par)
{
    A = Par[0];
	E = Par[1];
	Gamma = Par[2];
}

void Lorentzian::updatePar(double* C)
{
    A += AWF * C[0];
	E += EWF * C[1];
	Gamma += GWF * C[2];
}
