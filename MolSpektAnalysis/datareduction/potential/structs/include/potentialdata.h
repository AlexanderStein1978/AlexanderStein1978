//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POTENTIALDATA_H
#define POTENTIALDATA_H


#include "fitresult.h"
#include "constants.h"

struct SplinePoint;


struct PotentialData
{
	PotentialData();
	~PotentialData();
	
	int nPotCoeff, nLRCoeff, *pLRCoeff, NAdCorr, TAdCorr, PAdCorr, NSpinRGamma;
	double Ri, Ra, iExp, iOff, iCoeff, b, Rm, Tm, De, *PotCoeff, *LRCoeff, A_ex, beta, gamma, q_e, q_f, *adCorr;
	double RIso1AdCorr, RIso2AdCorr, AdCorr_b, AdCorr_Rm, *SpinRGamma, SpinRR1, SpinRR2;
	bool *CFree, *LRCFree, EIFree, *adCorrFree, shortRVariable;
	double *A, *C, Uinf,  iA, iO, alpha, *LRC, Salpha, Sbeta, *SA, *SLRC, Sinf, B, alpha1;
	int NA, N, NLRC, *PLRC;
	int p, q, NCi, Nb, *pCi, *pLRC, NSA, NSLRC, *pSLRC;
	double *bA, *Ci, Re, Rref, T, Ro, FQS;
	SplinePoint *points;
	PotentialType PotType;
	FitResult fitResult;
};

#endif
