//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "MLRPot.h"
#include "potentialdata.h"

#include <QMutex>
#include <QFile>
#include <QTextStream>


MLRPot::MLRPot(PotFit *Fit) : PotWorker(Fit, MorseLongRangePotential)
{
    Ci = beta = dCiLR = 0;
    pLRC = 0;
    p = q = NLRC = Nbeta = 0;
    De = Re = Rep = Rref = Rrefp = Rrefq = URe = betainf = TAs = dULRdRe = 0.0;
}

MLRPot::MLRPot(const MLRPot& C) : PotWorker(C)
{
    int n;
    if (C.NLRC > 0)
    {
        Ci = new double[NLRC = C.NLRC];
        pLRC = new int[NLRC];
        dCiLR = new double[NLRC];
        for (n=0; n < NLRC; n++) 
        {
            Ci[n] = C.Ci[n];
            pLRC[n] = C.pLRC[n];
            dCiLR[n] = C.dCiLR[n];
        }
    }
    else
    {
        Ci = 0;
        pLRC = 0;
        NLRC = 0;
    }
    if (C.Nbeta > 0)
    {
        beta = new double[Nbeta = C.Nbeta];
        for (n=0; n < Nbeta; n++) beta[n] = C.beta[n];
    }
    else
    {
        beta = 0;
        Nbeta = 0;
    }
    p = C.p;
    q = C.q;
    TAs = C.TAs;
    De = C.De;
    Re = C.Re;
    Rep = C.Rep;
    Rref = C.Rref;
    Rrefp = C.Rrefp;
    Rrefq = C.Rrefq;
    URe = C.URe;
    betainf = C.betainf;
    dULRdRe = C.dULRdRe;
}

MLRPot::~MLRPot()
{
    if (Ci != 0) 
    {
        delete[] Ci;
        delete[] pLRC;
        delete[] dCiLR;
    }
    if (beta != 0) delete[] beta;
}

void MLRPot::getCoefficients(int& rp, int& rq, double& rTAs, double& rDe, double& rRe, 
                             double& rRref, int& rNbeta, double*& rbeta, int& rNLRC, 
                             int *&rpLRC, double*& rLRC)
{
    int n;
    if (NLRC > 0)
    {
        rLRC = new double[rNLRC = NLRC];
        rpLRC = new int[NLRC];
        for (n=0; n < NLRC; n++) 
        {
            rLRC[n] = Ci[n];
            rpLRC[n] = pLRC[n];
        }
    }
    else
    {
        rLRC = 0;
        rpLRC = 0;
        rNLRC = 0;
    }
    if (Nbeta > 0)
    {
        rbeta = new double[rNbeta = Nbeta];
        for (n=0; n < Nbeta; n++) rbeta[n] = beta[n];
    }
    else
    {
        rbeta = 0;
        rNbeta = 0;
    }
    rp = p;
    rq = q;
    rTAs = TAs + De;
    rDe = De;
    rRe = Re;
    rRref = Rref;
}

void MLRPot::getDerivatives(double r, double& E, double* dbeta, double* dC, double& dRe, double& dDe)
{
    int n;
    double Rq = pow(r, q), zDe = 2.0 * De, UExp, myp;
    double Rp = Rq * pow(r, p-q), Rs;
    double yqRef = (Rq - Rrefq) / (Rq + Rrefq), ypRef = (Rp - Rrefp) / (Rp + Rrefp), yep = (Rp - Rep) / (Rs = Rp + Rep);
    double ULR, b = 0.0, mypRef = (1.0 - ypRef), ULRRe, Exp; 
    for (n=1, ULR = Ci[0] * (dC[0] = pow(r, -pLRC[0])); n < NLRC; n++)
    {
        dC[n] = dC[n-1] * pow(r, pLRC[n-1] - pLRC[n]);
        ULR += Ci[n] * dC[n];
    }
    for (n = Nbeta - 1; n>=0; n--)
    {
        b *= yqRef;
        b += beta[n];
    }
    b *= mypRef;
    b += ypRef * betainf;
    double zUM, mUExp = 1.0 - (UExp = (ULRRe = ULR * URe) * (Exp = exp(-b * yep)));
    E = TAs + De * mUExp * mUExp;
    double db = (zUM = zDe * UExp * mUExp) * yep * mypRef;
    for (n=1, dbeta[0] = db; n < Nbeta; n++) dbeta[n] = yqRef * dbeta[n-1];
    double pExp = zDe * mUExp * Exp * URe, pyepRef = ULRRe * (1.0 - (myp = yep * ypRef));
    for (n=0; n < NLRC; n++) dC[n] = pExp * (dCiLR[n] * pyepRef - dC[n]);
    dRe = zUM * (dULRdRe * (myp - 1.0) - b * double(p) * Rep / (Re * Rs) * (1.0 + yep)); 
    dDe = mUExp * mUExp + zUM * myp / De;
}

void MLRPot::getLRCoeff(double& R, int& numCoefficients, int*& Exponents, double*& Coefficients)
{
    R = 0.0;
    numCoefficients = NLRC;
    Exponents = pLRC;
    Coefficients = Ci;
}

void MLRPot::getNC(int& N, int& NC, int& nLRC)
{
    N = (NC = Nbeta) + (nLRC = NLRC) + 3; 
}

double MLRPot::GetDe()
{
    return De;
}

void MLRPot::getMinimum(double& R, double& T)
{
    Lock->lock();
    R = Re;
    T = TAs - De;
    Lock->unlock();
}

PotentialData* MLRPot::getPotentialData()
{
    int n;
    PotentialData *R = new PotentialData();
    Lock->lock();
    R->p = p;
    R->q = q;
    R->NCi = NLRC;
    R->Nb = Nbeta;
    R->pCi = new int[NLRC];
    R->NAdCorr = NAdCorr;
    R->TAdCorr = TAdCorr;
    R->PAdCorr = PAdCorr;
    R->bA = new double[Nbeta];
    for (n=0; n < Nbeta; n++) R->bA[n] = beta[n];
    R->Ci = new double[NLRC];
    for (n=0; n < NLRC; n++)
    {
        R->Ci[n] = LRC[n];
        R->pCi[n] = pLRC[n];
    }
    R->Re = Re;
    R->Rref = Rref;
    R->T = TAs;
    R->De = De;
    if (NAdCorr > 0)
    {
        R->adCorr = new double[NAdCorr];
        for (n=0; n < NAdCorr; n++) R->adCorr[n] = adCorr[n];
    }
    R->RIso1AdCorr = RIso1AdCorr;
    R->RIso2AdCorr = RIso2AdCorr;
    R->AdCorr_b = AdCorr_b;
    R->AdCorr_Rm = AdCorr_Rm;
    Lock->unlock();
    return R;
}

double MLRPot::Point(double r, int)
{
    int n;
    double Rq = pow(r, q);
    double Rp = Rq * pow(r, p-q);
    double yqRef = (Rq - Rrefq) / (Rq + Rrefq), ypRef = (Rp - Rrefp) / (Rp + Rrefp);
    double ULR = 0.0, b = 0.0;
    for (n = NLRC - 1; n>=0; n--)
    {
        ULR += Ci[n];
        ULR *= pow(r, (n>0 ? pLRC[n-1] - pLRC[n] : -pLRC[n]));
    }
    for (n = Nbeta - 1; n>=0; n--)
    {
        b *= yqRef;
        b += beta[n];
    }
    b *= (1.0 - ypRef);
    b += ypRef * betainf;
    b = 1.0 - ULR * URe * exp(- b * (Rp - Rep) / (Rp + Rep));
    //printf("exp=%e, ", pr);
    return TAs + De * b * b;
}

double MLRPot::getUinf()
{
    return TAs;
}

MLRPot* MLRPot::scalePotential(double newRe, double newDe)
{
    Lock->lock();
    MLRPot *rPot = new MLRPot(0);
    int n, *rpLRC = new int[NLRC];
    double RF = newRe / Re, DF = newDe / De;
    double *rbeta = new double[Nbeta], *rCi = new double[NLRC];
    for (n=0; n < Nbeta; n++) rbeta[n] = beta[n];
    for (n=0; n < NLRC; n++) rCi[n] = Ci[n] * DF * pow(RF, -(rpLRC[n] = pLRC[n]));
    rPot->setCoefficients(p, q, TAs * DF, newDe, newRe, Rref * RF, Nbeta, rbeta, NLRC, 
                          rpLRC, rCi);
    Lock->unlock();
    return rPot;
}

void MLRPot::setCoefficients(int np, int nq, double nTAs, double nDe, double nRe, double nRref,
                             int nNbeta, double* nbeta, int nNLRC, int* npLRC, double* nLRC)
{
    p = np;
    q = nq;
    TAs = nTAs - nDe;
    De = nDe;
    Re = nRe;
    Rref = nRref;
    Nbeta = nNbeta;
    if (beta != 0) delete[] beta;
    beta =  nbeta;
    if (nNLRC != NLRC) 
    {
      delete[] dCiLR;
      dCiLR = new double[nNLRC];
    }
    NLRC = nNLRC;
    if (pLRC != 0) delete[] pLRC;
    if (Ci != 0) delete[] Ci;
    Ci = nLRC;
    pLRC = npLRC;
    calcRIndP();
}

void MLRPot::calcRIndP()
{
    int n;
    for (n=0, URe = dULRdRe = 0.0; n < NLRC; n++) 
    {
        URe += Ci[n] * (dCiLR[n] = pow(Re, -pLRC[n]));
        dULRdRe += pLRC[n] * Ci[n] * dCiLR[n] / Re;
    }
    dULRdRe *= (URe = 1.0 / URe);
    Rep = pow(Re, p);
    Rrefp = pow(Rref, p);
    Rrefq = pow(Rref, q);
    betainf = log(2.0 * De * URe);
}

void MLRPot::shift(double Energy)
{
    TAs += Energy;
}

void MLRPot::test()
{
    QFile F("test.dat");
    F.open(QIODevice::WriteOnly);
    QTextStream S(&F);
    double B, V, R, E, dBeta[Nbeta], dC[NLRC], dRe, dDe, ds, lE;
    int n;
    for (R = 0.5 * Re; R <= 5 * Re; R += (R < Re ? 0.1 * Re : Re))
    {
        S << "R = " << QString::number(R, 'f', 2) << ":\n";
        for (n=0; n < Nbeta; n++)
        {
            S << "beta[" << QString::number(n) << "]:\n";
            ds = 100.0 / (B = beta[n]);
            for (V = 0.9 * B, lE = 0.0; (V > 0.0 ? V <= 1.1 * B : V >= 1.1 * B); V += 0.01 * B, lE = E)
            {
                beta[n] = V;
                calcRIndP();
                getDerivatives(R, E, dBeta, dC, dRe, dDe);
                S << "V=" << QString::number(V, 'g', 3) << ", E=" << QString::number(E, 'g', 7)
                  << ", dE/dbeta=" << QString::number(dBeta[n], 'g', 7) << " ~ " 
                  << (lE != 0.0 ? QString::number((E - lE) * ds, 'g', 7) : "?") << "\n";
            }
            beta[n] = B;
        }
        for (n=0; n < NLRC; n++)
        {
            S << "Ci[" << QString::number(n) << "]:\n";
            ds = 100 / Ci[n];
            for (V = 0.9 * (B = Ci[n]), lE = 0.0; (V > 0.0 ? V <= 1.1 * B : V >= 1.1 * B); V += 0.01 * B, lE = E)
            {
                Ci[n] = V;
                calcRIndP();
                getDerivatives(R, E, dBeta, dC, dRe, dDe);
                S << "V=" << QString::number(V, 'g', 3) << ", E=" << QString::number(E, 'g', 7)
                  << ", dE/dC=" << QString::number(dC[n], 'g', 7) << " ~ " 
                  << (lE != 0.0 ? QString::number((E - lE) * ds, 'g', 7) : "?") << "\n";
            }
            Ci[n] = B;
        }
        S << "Re:\n";
        ds = 100 / Re;
        for (V = 0.9 * (B = Re), lE = 0.0; V <= 1.1 * B; V += 0.01 * B, lE = E)
        {
            Re = V;
            calcRIndP();
            getDerivatives(R, E, dBeta, dC, dRe, dDe);
            S << "V=" << QString::number(V, 'g', 3) << ", E=" << QString::number(E, 'g', 7)
              << ", dE/dRe=" << QString::number(dRe, 'g', 7) << " ~ " 
              << (lE != 0.0 ? QString::number((E - lE) * ds, 'g', 7) : "?") << "\n";
        }
        Re = B;
        S << "De:\n";
        ds = 100 / De;
        for (V = 0.9 * (B = De), lE = 0.0; V <= 1.1 * B; V += 0.01 * B, lE = E)
        {
            De = V;
            calcRIndP();
            getDerivatives(R, E, dBeta, dC, dRe, dDe);
            S << "V=" << QString::number(V, 'g', 3) << ", E=" << QString::number(E, 'g', 7)
              << ", dE/dDe=" << QString::number(dDe, 'g', 7) << " ~ " 
              << (lE != 0.0 ? QString::number((E - lE) * ds, 'g', 7) : "?") << "\n";
        }
        De = B;
    }
    calcRIndP();
}
