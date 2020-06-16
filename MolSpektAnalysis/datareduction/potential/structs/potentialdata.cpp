//
// C++ Implementation: PotentialData
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "potentialdata.h"
#include "SplinePoint.h"


PotentialData::PotentialData()
{
    nPotCoeff = nLRCoeff = 0;
    pLRCoeff = 0;
    NAdCorr = TAdCorr = PAdCorr = NSpinRGamma = 0;
    Ri = Ra = Ro = iExp = iOff = iCoeff = b = Rm = Tm = De = 0.0;
    PotCoeff = LRCoeff = 0;
    A_ex = beta = gamma = q_e = q_f = 0.0;
    adCorr = 0;
    RIso1AdCorr = RIso2AdCorr = 0.0;
    SpinRGamma = 0;
    SpinRR1 = SpinRR2 = 0.0;
    CFree = LRCFree = 0;
    EIFree = true;
    adCorrFree = 0;
    shortRVariable = true;
    A = C = 0;
    Uinf = iA = iO = alpha = 0.0;
    LRC = 0;
    Salpha = Sbeta = 0.0;
    SA = SLRC = 0;
    Sinf = B = alpha1 = 0.0;
    NA = N = NLRC = 0;
    PLRC = 0;
    p = q = NCi = Nb = 0;
    pCi = 0;
    pLRC = 0;
    NSA = NSLRC = 0;
    pSLRC = 0;
    bA = Ci = 0;
    Re = Rref = T = De = 0.0;
    AdCorr_b = AdCorr_Rm = 0.0;
    points = 0;
    PotType = NoPotential;
}

PotentialData::~PotentialData()
{
    if (pLRCoeff != 0) delete[] pLRCoeff;
    if (PotCoeff != 0) delete[] PotCoeff;
    if (LRCoeff != 0) delete[] LRCoeff;
    if (adCorr != 0) delete[] adCorr;
    if (SpinRGamma != 0) delete[] SpinRGamma;
    if (CFree != 0) delete[] CFree;
    if (LRCFree != 0) delete[] LRCFree;
    if (adCorrFree != 0) delete[] adCorrFree;
    if (A != 0) delete[] A;
    if (C != 0) delete[] C;
    if (LRC != 0) delete[] LRC;
    if (SA != 0) delete[] SA;
    if (SLRC != 0) delete[] SLRC;
    if (PLRC != 0) delete[] PLRC;
    if (pCi != 0) delete[] pCi;
    if (pLRC != 0) delete[] pLRC;
    if (pSLRC != 0) delete[] pSLRC;
    if (bA != 0) delete[] bA;
    if (Ci != 0) delete[] Ci;
    if (points !=0 ) delete[] points;
}
