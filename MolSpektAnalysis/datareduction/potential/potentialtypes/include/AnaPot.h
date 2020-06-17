//
// C++ Interface: AnaPot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef ANAPOT_H
#define ANAPOT_H


#include "potworker.h"


class AnaPot : public PotWorker
{
	public:
		AnaPot(PotFit *Fit);
		AnaPot(const AnaPot &C);
		~AnaPot();
		AnaPot *scalePotential(double newRe, double newDe);
		void setInnerWall(double R, double Offset, double Exponent, double Coefficient);
		void getInnerWall(double &R, double &Offset, double &Exponent, double &Coefficient);
		void setPotCoeff(double b, double Rm, double Tm, double De, int numCoefficients, 
						 double *Coefficients, bool *CoeffFree = 0);
		void getPotCoeff(double &b, double &Rm, double &Tm, double &De, int &numCoefficients,
						 double *&Coefficients);
		void setLRCoeff(double R, int numCoefficients, int *Exponents, double *Coefficients, bool *LRCFree = 0);
		void setExchangeCoeff(double A_ex, double beta, double gamma, bool Free = true);
		void getExchangeCoeff(double &A_ex, double &beta, double &gamma);
		void cdConnectLR1C();
		void cdConnectLR(int firstCoeff);
		void cdConnectSR();
		double calcPotCoeff(int numPoints, double *R, double *U, double *Sig, int numCoeff, double Rm, double b);
		double fitPotential(int NPoints, double R0, double h, double *U, double *Sig, int NCoeff, int NLRC, int *pLRC,
							double Ri, double Ra, double Rm, double b, double Uinf, double iExp);
		void fitAdCorrToRmb();
		void getMinimum(double &R, double &T);
		double getUinf();
		void sConnect();
		void shift(double Energy);
		void getConR(double &Ri, double &Ra);
		void getNC(int &N, int &NC, int &NLRC);
		void getDerivatives(double r, double *dC, double *dLRC);
		PotentialData *getPotentialData();
		
		inline void setDe(double nDe)
		{
			Uinf = nDe;
		}
		
		inline bool *getCoeffFree()
		{
			return CoeffFree;
		}
		
		inline bool *getLRCFree()
		{
			return LRCFree;
		}
		
		inline bool isExcIntFree()
		{
			return ExcIntFree;
		}
		
		inline bool isExcIntAvailable()
		{
			return ExcIntFree || A_ex != 0.0;
		}
		
		inline void getRmb(double &rRm, double &rb)
		{
			rRm = Rm;
			rb = b;
		}
		
		inline double getIExp()
		{
			return iExp;
		}
		
	protected:
		double *Points(double Rmin, double Rmax, int nPoints, int FC = 0);
		void saveCoefficients(double *&bC);
		void updatePotential(double *C);
		void getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients);
		void createSMap();
		void restorePotential(double *bC);
		void calc_hCi_cP_EQS(int NEQ, double **EQS, double *Sig, double *EDiff, double hCi_cP);
		void calcS(int numPoints, double Rmin, double Rmax);
		double Point(double r, int typ = 0);
		int cAdRows();
		
	private:
		double Ri, Ra, iOff, iCoeff, *PotCoeff, *LRCoeff, Tm, *SWF;
		int nLRCoeff, *pLRCoeff;
		int nPotCoeff, numFreeCoeff, numFreeLRC;
		double Rm, b, A_ex, beta, gamma;
		bool *LRCFree, *CoeffFree, ExcIntFree, ADLRC;
};

#endif
