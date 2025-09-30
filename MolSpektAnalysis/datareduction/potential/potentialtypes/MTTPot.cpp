//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "MTTPot.h"
#include "potentialdata.h"

#include <QMutex>


MTTPot::MTTPot(PotFit *Fit) : PotWorker(Fit, ModifiedTangToenniesPotential)
{
    A = LRC = C = CLRC = 0;
    NA = NLRC = NC = NCLRC = 0;
    alpha = beta = gamma = Uinf = Calpha = Cbeta = CInf = B = alpha1 = 0.0;
    pLRC = pCLRC = 0;
}

MTTPot::MTTPot(const MTTPot& Co) : PotWorker(Co)
{
    int n;
    if ((NA = Co.NA) > 0) for (n=0, A = new double[NA]; n < NA; n++) A[n] = Co.A[n];
    else A = 0;
    if ((NLRC = Co.NLRC) > 0) 
        for (n=0, LRC = new double[NLRC], pLRC = new int[NLRC]; n < NLRC; n++)
    {
        LRC[n] = Co.LRC[n];
        pLRC[n] = Co.pLRC[n];
    }
    B = Co.B;
    alpha1 = Co.alpha1;
    alpha = Co.alpha;
    beta = Co.beta;
    gamma = Co.gamma;
    Uinf = Co.Uinf;
    if ((NC = Co.NC) > 0) for (n=0, C = new double[NC]; n < NC; n++) C[n] = Co.C[n];
    else C = 0;
    if ((NCLRC = Co.NCLRC) > 0) 
        for (n=0, CLRC = new double[NCLRC], pCLRC = new int[NCLRC]; n < NCLRC; n++)
    {
        CLRC[n] = Co.CLRC[n];
        pCLRC[n] = Co.pCLRC[n];
    }
    Calpha = Co.Calpha;
    Cbeta = Co.Cbeta;
    CInf = Co.CInf;
}

MTTPot::~MTTPot()
{
    if (NA > 0) delete[] A;
    if (NLRC > 0)
    {
        delete[] LRC;
        delete[] pLRC;
    }
    if (NC > 0) delete[] C;
    if (NCLRC > 0)
    {
        delete[] CLRC;
        delete[] pCLRC;
    }
}

void MTTPot::getCoefficients(double &rB, double &ralpha1, int& rNA, double*& rA, double& rAlpha, double& rBeta, 
                             double& rGamma, int& rNLRC, int*& rpLRC, double*& rLRC, 
                             double &rUinf, int &rNC, double *&rC, double &rCalpha, 
                             double &rCbeta, int &rNCLRC, int *&rpCLRC, double *&rCLRC, 
                             double &rCInf)
{
    int n;
    if ((rNA = NA) > 0) for (n=0, rA = new double[NA]; n < NA; n++) rA[n] = A[n];
    else rA = 0;
    rB = B;
    ralpha1 = alpha1;
    rAlpha = alpha;
    rBeta = beta;
    rGamma = gamma;
    if ((rNLRC = NLRC) > 0) 
        for (n=0, rLRC = new double[NLRC], rpLRC = new int[NLRC]; n < NLRC; n++)
    {
        rpLRC[n] = pLRC[n];
        rLRC[n] = LRC[n];
    }
    else 
    {
        rLRC = 0;
        rpLRC = 0;
    }
    rUinf = Uinf;
    if ((rNC = NC) > 0) for (n=0, rC = new double[NC]; n < NC; n++) rC[n] = C[n];
    else rC = 0;
    rCalpha = Calpha;
    rCbeta = Cbeta;
    if ((rNCLRC = NCLRC) > 0)
        for (n=0, rCLRC = new double[NCLRC], rpCLRC = new int[NCLRC]; n < NCLRC; n++)
    {
        rpCLRC[n] = pCLRC[n];
        rLRC[n] = CLRC[n];
    }
    rCInf = CInf;
}

double MTTPot::GetDe(int c)
{
    double R, V;
    getMinimum(R, V, c);
    switch (c)
    {
        case 0:
            return Uinf - V;
        case -1:
            return Uinf - V - CInf;
        case 1:
            return Uinf - V + 2.0 * CInf;
            break;
    }
    return 0.0;
}

void MTTPot::getMinimum(double& R, double& V, int c)
{
    double MD = 1e-5, CMin = Uinf - 1e-4, MinR = -1.0, LMin = Uinf, S, Ri = rmin, Ra = rmax;
    Lock->lock();
    for (S=1.0; LMin - CMin > MD; S *= 0.1)
    {
        for (R = Ri, LMin = CMin; R <= Ra; R+=S) if ((V = Point(R, c)) < CMin)
        {
            CMin = V;
            MinR = R;
        }
        Ri = R - S;
        Ra = R + S;
    }
    R = MinR;
    V = CMin;
    Lock->unlock();
}

PotentialData* MTTPot::getPotentialData()
{
    int n;
    PotentialData *R = new PotentialData();
    Lock->lock();
    R->NA = NA;
    R->NLRC = NLRC;
    if (NLRC > 0)
    {
        R->pLRC = new int[NLRC];
        R->LRC = new double[NLRC];
        for (n=0; n < NLRC; n++)
        {
            R->pLRC[n] = pLRC[n];
            R->LRC[n] = LRC[n];
        }
    }
    R->NSA = NC;
    R->NSLRC = NCLRC;
    R->alpha = alpha;
    R->beta = beta;
    R->gamma = gamma;
    if (NA > 0)
    {
        R->A = new double[NA];
        for (n=0; n < NA; n++) R->A[n] = A[n];
    }
    R->Uinf = Uinf;
    R->Salpha = Calpha;
    R->Sbeta = Cbeta;
    if (NC > 0)
    {
        R->SA = new double[NC];
        for (n=0; n < NC; n++) R->SA[n] = C[n];
    }
    if (NCLRC > 0)
    {
        R->pSLRC = new int[NCLRC];
        R->SLRC = new double[NCLRC];
        for (n=0; n < NCLRC; n++)
        {
            R->pSLRC[n] = pCLRC[n];
            R->SLRC[n] = CLRC[n];
        }
    }
    R->Sinf = CInf;
    R->B = B;
    R->alpha1 = alpha1;
    Lock->unlock();
    return R;
}

double MTTPot::Point(double r, int c)
{
    double bR = beta * r, V = 0.0, FS = 1.0, dR = 1.0 / r, RP = 1.0, F = 1.0, BRP = 1.0;
    double expBR = exp(-bR);
    int n, m, m1, m2 = 0;
    for (n = NA - 1; n >= 0; n--)
    {
        V *= r;
        V += A[n];
    }
    V *= exp((-alpha - gamma * r) * r);
    V += Uinf + B * exp(-alpha1 * r);
    for (n=0; n < NLRC; n++)
    {
        m1 = m2;
        m2 = pLRC[n];
        for (m = m1 + 1; m <= m2; m++)
        {
            BRP *= bR;
            F *= m;
            FS += BRP / F;
            RP *= dR;
        }
        V -= LRC[n] * RP * (1.0 - FS * expBR);
    }
    if (c != 0)
    {
        bR = Cbeta * r, FS = 0.0, RP = 1.0, F = 1.0, BRP = 1.0;
        expBR = exp(-bR);
        for (n = NC - 1; n >= 0; n--)
        {
            FS *= r;
            FS += C[n];
        }
        FS *= exp(-Calpha * r);
        if (c==-1) V -= FS + CInf;
        else V += 2.0 * (FS + CInf);
        for (n = m2 = 0, FS = 1.0; n < NCLRC; n++)
        {
            m1 = m2;
            m2 = pCLRC[n];
            for (m = m1 + 1; m <= m2; m++)
            {
                BRP *= bR;
                F *= m;
                FS += BRP / F;
                RP *= dR;
            }
            if (c==-1) V += CLRC[n] * RP * (1.0 - FS * expBR);
            else V -= 2.0 * CLRC[n] * RP * (1.0 - FS * expBR);
        }
    }
    return V;
}

void MTTPot::getLRCoeff(double& R, int& numCoefficients, int*& Exponents, double*& Coefficients)
{
    R = 0.0;
    numCoefficients = NC;
    Exponents = pCLRC;
    Coefficients = C;
}

double MTTPot::GetUinf(int c)
{
    switch (c)
    {
        case 0:
            return Uinf;
        case -1:
            return Uinf - CInf;
        case 1:
            return Uinf + 2.0 * CInf;
            break;
    }
    return 0.0;
}

MTTPot* MTTPot::scalePotential(double newRe, double newDe)
{
    double Re, De;
    int n;
    getMinimum(Re, De);
    Lock->lock();
    De = Uinf - De;
    double XSF = Re / newRe, YSF = newDe / De;
    MTTPot *RPot = new MTTPot(0);
    RPot->alpha = XSF * alpha;
    RPot->alpha1 = XSF * alpha1;
    RPot->B = YSF * B;
    RPot->gamma = XSF * XSF * gamma;
    RPot->A = new double[RPot->NA = NA];
    for (n=0; n < NA; n++) RPot->A[n] = YSF * pow(XSF, n) * A[n];
    RPot->LRC = new double[RPot->NLRC = NLRC];
    RPot->pLRC = new int[NLRC];
    for (n=0; n < NLRC; n++)
    {
        RPot->pLRC[n] = pLRC[n];
        RPot->LRC[n] = YSF * pow(XSF, -pLRC[n]) * LRC[n];
    }
    RPot->beta = XSF * beta;
    RPot->Uinf = Uinf;
    RPot->Calpha = XSF * Calpha;
    if ((RPot->NC = NC) > 0)
    {
        RPot->C = new double[NC];
        for (n=0; n < NC; n++) RPot->C[n] = YSF * pow(XSF, n) * C[n];
    }
    else RPot->C = 0;
    if ((RPot->NCLRC = NCLRC) > 0)
    {
        RPot->CLRC = new double[NCLRC];
        RPot->pCLRC = new int[NCLRC];
        for (n=0; n < NCLRC; n++)
        {
            RPot->pCLRC[n] = pCLRC[n];
            RPot->CLRC[n] = YSF * pow(XSF, -pCLRC[n]) * CLRC[n];
        }
    }
    else
    {
        RPot->pCLRC = 0;
        RPot->CLRC = 0;
    }
    RPot->Cbeta = XSF * beta;
    RPot->CInf = CInf;
    Lock->unlock();
    return RPot;
}

void MTTPot::setCoefficients(double nB, double nalpha1, int nNA, double* nA, double nAlpha, double nBeta, 
                             double nGamma, int nNLRC, int* npLRC, double* nLRC, 
                             double nUinf, int nNC, double *nC, double nCalpha, 
                             double nCbeta, int nNCLRC, int *npCLRC, double *nCLRC, 
                             double nCInf)
{
    if (NA > 0) delete[] A;
    if (NLRC > 0)
    {
        delete[] LRC;
        delete[] pLRC;
    }
    if (NC > 0) delete[] C;
    if (NCLRC > 0)
    {
        delete[] CLRC;
        delete[] pCLRC;
    }
    B = nB;
    alpha1 = nalpha1;
    NA = nNA;
    A = nA;
    alpha = nAlpha;
    beta = nBeta;
    gamma = nGamma;
    NLRC = nNLRC;
    pLRC = npLRC;
    LRC = nLRC;
    Uinf = nUinf;
    NC = nNC;
    C = nC;
    Calpha = nCalpha;
    Cbeta = nCbeta;
    NCLRC = nNCLRC;
    pCLRC = npCLRC;
    CLRC = nCLRC;
    CInf = nCInf;
}
