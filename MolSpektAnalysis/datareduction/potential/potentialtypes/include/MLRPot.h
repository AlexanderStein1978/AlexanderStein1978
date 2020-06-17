//
// C++ Interface: MLRPot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MLRPOT_H
#define MLRPOT_H


#include "potworker.h"


class MLRPot : public PotWorker
{
	public:
		MLRPot(PotFit *Fit);
		MLRPot(const MLRPot &C);
		~MLRPot();
		MLRPot *scalePotential(double newRe, double newDe);
		void getMinimum(double &R, double &T);
		void getDerivatives(double R, double &E, double *dbeta, double *dC, double &dRe, double &dDe);
		void getNC(int &N, int &NC, int &NLRC);
		double GetDe();
		double getUinf();
		void shift(double Energy);
		void setCoefficients(int p, int q, double TAs, double De, double Re, double Rref, 
							 int Nbeta, double *beta, int NLRC, int *pLRC, double *LRC);
		void getCoefficients(int &p, int &q, double &TAs, double &De, double &Re, double &Rref,
							 int &Nbeta, double *&beta, int &NLRC, int *&pLRC, double *&LRC);
		void calcRIndP();
		void test();
		PotentialData *getPotentialData();
		
	protected:
		void getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients);
		double Point(double r, int FC = 0);
		
	private:
		int Nbeta, NLRC, p, q, *pLRC;
		double Re, Rep, Rref, Rrefq, Rrefp, De, *Ci, *beta, betainf, URe, TAs, *dCiLR, dULRdRe;
};

#endif
