//
// C++ Implementation: AnaPot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "AnaPot.h"
#include "tools.h"
#include "utils.h"
#include "potentialdata.h"

#include <QMutex>


AnaPot::AnaPot(PotFit *Fit) : PotWorker(Fit, analyticalPotential)
{
    nLRCoeff = nPotCoeff = numFreeCoeff = numFreeLRC = 0;
    b = iCoeff = iOff = Ra = Ri = iExp = A_ex = beta = gamma = Rm = Tm = 0.0;
    LRCoeff = PotCoeff = 0;
    pLRCoeff = 0;
    SWF = 0;
    LRCFree = CoeffFree = 0;
    ExcIntFree = false;
    ADLRC = false;
}

AnaPot::AnaPot(const AnaPot &C) : PotWorker(C)
{
    //printf("AnaPot::AnaPot(const AnaPot &C)\n");
    int n;
    if (C.SWF != 0)
    {
        SWF = new double[C.nPotCoeff + C.nLRCoeff - 1];
        for (n=0; n < C.nPotCoeff + C.nLRCoeff - 1; n++) SWF[n] = C.SWF[n];
    }
    else SWF = 0;
    b = C.b;
    iExp = C.iExp;
    numFreeCoeff = C.numFreeCoeff;
    numFreeLRC = C.numFreeLRC;
    ExcIntFree = C.ExcIntFree;
    if ((nLRCoeff = C.nLRCoeff) > 0)
    {
        LRCoeff = new double[nLRCoeff + 2];
        pLRCoeff = new int[nLRCoeff + 2];
        LRCFree = new bool[nLRCoeff];
        for (n=0; n < nLRCoeff; n++) 
        {
            LRCoeff[n] = C.LRCoeff[n];
            pLRCoeff[n] = C.pLRCoeff[n];
            LRCFree[n] = C.LRCFree[n];
        }
    }
    else
    {
        LRCoeff = 0;
        pLRCoeff = 0;
        LRCFree = 0;
    }
    if ((nPotCoeff = C.nPotCoeff) > 0)
    {
        PotCoeff = new double[nPotCoeff];
        CoeffFree = new bool[nPotCoeff];
        for (n=0; n < nPotCoeff; n++) 
        {
            PotCoeff[n] = C.PotCoeff[n];
            CoeffFree[n] = C.CoeffFree[n];
        }
    }
    iCoeff = C.iCoeff;
    iOff = C.iOff;
    Ra = C.Ra;
    Ri = C.Ri;
    Rm = C.Rm;
    Tm = C.Tm;
    A_ex = C.A_ex;
    beta = C.beta;
    gamma = C.gamma;
    ADLRC = C.ADLRC;
    //printf("Ende AnaPot\n");
}

AnaPot::~AnaPot()
{
    if (nPotCoeff > 0) 
    {
        delete[] PotCoeff;
        delete[] CoeffFree;
    }
    if (nLRCoeff > 0)
    {
        delete[] pLRCoeff;
        delete[] LRCoeff;
        delete[] LRCFree;
    }
    if (SWF != 0) delete[] SWF;
}

AnaPot* AnaPot::scalePotential(double newRe, double newDe)
{
    AnaPot *RPot = new AnaPot(Fit);
    int n;
    double ARe, ADe, xSF, ySF;
    getMinimum(ARe, ADe);
    xSF = newRe / ARe;
    ySF = newDe / (Uinf - ADe);
    RPot->b = b;
    RPot->Uinf = ySF * Uinf;
    if ((RPot->nLRCoeff = nLRCoeff) > 0)
    {
        RPot->LRCoeff = new double[nLRCoeff + 2];
        RPot->pLRCoeff = new int[nLRCoeff + 2];
        for (n=0; n < nLRCoeff; n++) 
            RPot->LRCoeff[n] = LRCoeff[n] * pow(xSF, RPot->pLRCoeff[n] = pLRCoeff[n]) * ySF;
    }
    else
    {
        RPot->LRCoeff = 0;
        RPot->pLRCoeff = 0;
    }
    if ((RPot->nPotCoeff = nPotCoeff) > 0)
    {
        RPot->PotCoeff = new double[nPotCoeff];
        for (n=0; n < nPotCoeff; n++) RPot->PotCoeff[n] = ySF * PotCoeff[n];
    }
    RPot->iCoeff = iCoeff * pow(xSF, RPot->iExp = iExp) * ySF;
    RPot->iOff = ySF * iOff;
    RPot->Ra = xSF * Ra;
    RPot->Ri = xSF * Ri;
    RPot->Rm = xSF * Rm;
    RPot->Tm = ySF * Tm;
    RPot->A_ex = ySF * A_ex * pow(xSF, -gamma);
    RPot->beta = beta / xSF;
    RPot->gamma = gamma;
    return RPot;
}

double *AnaPot::Points(double RMin, double RMax, int nP, int, std::function<double(const double)> mapping)
{
    //printf("getPointRep: nP=%d\n", nP);
    if (nLRCoeff <= 0 || nPotCoeff <= 0 || Ri == 0.0 || nP <= 0 || RMin >= RMax) return NULL;
    double *RA = new double[nP], h = (RMax - RMin) / double(nP - 1), R = RMin, x, edR;
    int i, n, p;
    //printf("Ni=%d, Na=%d, Np=%d\n", Ni, Na, nP);
    if (mapping(RMin) > mapping(RMax))
    {
        R = RMax;
        h *= -1.0;
    }
    for (i=0; i < nP; i++)
    {
        const double r = mapping(R);
        if (r > Ri) break;
        RA[i] = iOff + iCoeff * pow(r, -iExp);
        //printf("RA[%d]=%g\n", i, RA[i]);
        R += h;
    }
    for (; i < nP; i++)
    {
        const double r = mapping(R);
        if (r > Ra) break;
        x = (r - Rm) / (r + b * Rm);
        for (n = nPotCoeff - 1, RA[i] = 0.0; n >= 0; n--)
        {
            RA[i] += PotCoeff[n];
            RA[i] *= x;
        }
        RA[i] += Tm;
        R += h;
    }
    for (; i < nP; i++)
    {
        const double r = mapping(R);
        edR = 1.0 / r;
        for (n = nLRCoeff - 1, p = pLRCoeff[nLRCoeff-1], RA[i] = 0.0; n >= 0; n--)
        {
            while (p > pLRCoeff[n]) 
            {
                RA[i] *= edR;
                p--;
            }
            RA[i] -= LRCoeff[n];
        }
        RA[i] *= pow(edR, pLRCoeff[0]);
        RA[i] += Uinf;
        if (A_ex != 0.0) RA[i] -= A_ex * pow(r, gamma) * exp(-beta * r);
        R += h;
    }
    if (h < 0.0) for (i=0; i < nP / 2; ++i)
    {
        const double B = RA[i];
        RA[i] = RA[nP - i - 1];
        RA[nP - i - 1] = B;
    }
    //printf("Ende getPointRep\n");
    return RA;
}

double AnaPot::Point(double r, int typ)
{
    if (typ == 0 ? r <= Ri : typ == 1) return iOff + iCoeff * pow(r, -iExp);
    int n, p;
    double R, x;
    if (typ == 0 ? r < Ra : typ == 2)
    {
        if (nPotCoeff <= 0) return 0.0;
        x = (r - Rm) / (r + b * Rm);
        for (n = nPotCoeff - 1, R = 0.0; n >= 0; n--)
        {
            R += PotCoeff[n];
            R *= x;
        }
        return R + Tm;
    }
    if (nLRCoeff <= 0) return 0.0;
    double edR = 1.0 / r;
    for (n = nLRCoeff - 1, p = pLRCoeff[nLRCoeff-1], R = 0.0; n >= 0; n--)
    {
        while (p > pLRCoeff[n])
        {
            R *= edR;
            p--;
        }
        R -= LRCoeff[n];
    }
    if (A_ex == 0.0) return R * pow(edR, pLRCoeff[0]) + Uinf;
    return R * pow(edR, pLRCoeff[0]) + Uinf - A_ex * pow(r, gamma) * exp(-beta * r);
}

void AnaPot::getMinimum(double &R, double &T)
{
    //printf("AnaPot::GetMin()\n");
    Lock->lock();
    double x1, y1, x2, y2=0.0, x3 = rmin, y3=Point(rmin), r;
    for (r = rmin + 0.1, T = y3; r <= rmax; r+=0.1)
    {
        x1 = x2;
        y1 = y2;
        x2 = x3;
        y2 = y3;
        x3 = r;
        y3 = Point(r);
        if (y1 != 0.0 && y1 >= y2 && y2 < y3 && y2 < T)
        {
            R = x2;
            T = y2;
        }
    }
    x2 = R;
    y2 = T;
    x1 = R - 0.1;
    y1 = Point(x1);
    x3 = R + 0.1;
    y3 = Point(x3);
    while (fabs(r - T) > 0.0001)
    {
        ParabInterpol(x1, y1, x2, y2, x3, y3, R, T);
        r = Point(R);
        if (R < x2)
        {
            x3 = x2;
            y3 = y2;
        }
        else
        {
            x1 = x2;
            y1 = y2;
        }
        x2 = R;
        y2 = r;
    }
    T = r;
    //printf("Ende GetMin, T=%f\n", T);
    Lock->unlock();
}

void AnaPot::setInnerWall(double R, double O, double E, double C)
{
    Ri = R;
    iOff = O;
    iExp = E;
    iCoeff = C;
}

void AnaPot::getInnerWall(double &R, double &O, double &E, double &C)
{
    R = Ri;
    O = iOff;
    E = iExp;
    C = iCoeff;
}

void AnaPot::getConR(double &i, double &a)
{
    i = Ri;
    a = Ra;
}

void AnaPot::getDerivatives(double, double*, double*)
{

}

void AnaPot::setLRCoeff(double R, int N, int *E, double *C, bool *lRCFree)
{
    int n, p=1;
    double Co;
    bool Free;
    numFreeCoeff -= numFreeLRC;
    if (nLRCoeff > 0)
    {
        if (E != pLRCoeff) delete[] pLRCoeff;
        if (C != LRCoeff) delete[] LRCoeff;
        if (lRCFree != LRCFree && (lRCFree != 0 || N != nLRCoeff))
        {
            delete[] LRCFree;
            LRCFree = 0;
        }
    }
    Ra = R;
    nLRCoeff = N;
    pLRCoeff = E;
    LRCoeff = C;
    if (lRCFree != 0) LRCFree = lRCFree;
    else if (LRCFree == 0 && N>0) for (n=0, LRCFree = new bool[N]; n<N; n++) LRCFree[n] = true;
    while (p!=0) for (n=1, p=0; n<N; n++) if (pLRCoeff[n] < pLRCoeff[n-1])
    {
        p = pLRCoeff[n];
        pLRCoeff[n] = pLRCoeff[n-1];
        pLRCoeff[n-1] = p;
        Co = LRCoeff[n-1];
        LRCoeff[n-1] = LRCoeff[n];
        LRCoeff[n] = Co;
        Free = LRCFree[n];
        LRCFree[n] = LRCFree[n-1];
        LRCFree[n-1] = Free;
    }
    for (numFreeLRC = n = 0; n < N-1; n++) if (LRCFree[n]) numFreeLRC++;
    numFreeCoeff += numFreeLRC;
    ADLRC = (N>4 ? E[N-1] == 16 && E[N-4] == 10 : false);
}

void AnaPot::getLRCoeff(double &R, int &N, int *&E, double *&C)
{
    R = Ra;
    N = nLRCoeff;
    E = pLRCoeff;
    C = LRCoeff;
}

void AnaPot::getNC(int& N, int& NC, int& NLRC)
{
    N = (NC = nPotCoeff) + (NLRC = nLRCoeff);
}

void AnaPot::setExchangeCoeff(double A, double b, double g, bool Free)
{
    A_ex = A;
    beta = b;
    gamma = g;
    ExcIntFree = Free;
}

void AnaPot::getExchangeCoeff(double& A, double& b, double& g)
{
    A = A_ex;
    b = beta;
    g = gamma;
}

void AnaPot::setPotCoeff(double nb, double nRm, double nTm, double nDe, int N, double *C, bool *Free)
{
    //printf("AnaPot::setPotCoeff: N=%d\n", N);
    int n;
    if (Free != 0) 
    {
        for (n = numFreeCoeff = 0; n<N; n++) if (Free[n]) numFreeCoeff++;
    }
    else if (N != nPotCoeff) numFreeCoeff = N;
    numFreeCoeff += numFreeLRC;
    if (nPotCoeff > 0)
    {
        if (C != PotCoeff) delete[] PotCoeff;
        if (Free != CoeffFree && (Free != 0 || N != nPotCoeff))
        {
            delete[] CoeffFree;
            CoeffFree = 0;
        }
    }
    b = nb;
    Rm = nRm;
    Tm = nTm;
    Uinf = nDe;
    nPotCoeff = N;
    PotCoeff = C;
    if (Free != 0) CoeffFree = Free;
    else if (CoeffFree == 0 && N>0) for (CoeffFree = new bool[N], n=0; n<N; n++) CoeffFree[n] = true;
}

void AnaPot::getPotCoeff(double &nb, double &nRm, double &nTm, double &nDe, int &N, double *&C)
{
    nb = b;
    nRm = Rm;
    nTm = Tm;
    nDe = Uinf;
    N = nPotCoeff;
    C = PotCoeff;
    //printf("AnaPot::getPotCoeff: N=%d\n", N);
}

void AnaPot::cdConnectLR1C()
{
    int n, p;
    double dxdR, z, A = 0.0, dAdR = 0.0, x, L = 0.0, dLdR = 0.0, edR, D, pnC;
    /*if (ADLRC) for (nLRCoeff = 1; pLRCoeff[nLRCoeff - 1] < 10; nLRCoeff++) ;
    else if (nLRCoeff >= 3 ? pLRCoeff[nLRCoeff - 3] == 6 && pLRCoeff[nLRCoeff - 2] == 8 
                && pLRCoeff[nLRCoeff - 1] == 10 : false)
    {
        ADLRC = true;
        int *pB = new int[nLRCoeff + 3];
        double *LB = new double[nLRCoeff + 3];
        bool *fB = new bool[nLRCoeff + 3];
        for (n=0; n < nLRCoeff; n++)
        {
            pB[n] = pLRCoeff[n];
            LB[n] = LRCoeff[n];
            fB[n] = LRCFree[n];
        }
        delete[] pLRCoeff;
        delete[] LRCoeff;
        delete[] LRCFree;
        pLRCoeff = pB;
        LRCoeff = LB;
        LRCFree = fB;
    }*/
    Lock->lock();
    z = Ra + b * Rm;
    x = (Ra - Rm) / z;
    dxdR = (1.0 - x) / z;
    for (n = nPotCoeff - 1; n >= 0; n--)
    {
        dAdR += double(n+1) * PotCoeff[n];
        if (n>0) dAdR *= x;
        A += PotCoeff[n];
        A *= x;
    }
    dAdR *= dxdR;
    A += Tm;
    edR = 1.0 / Ra;
    if (nLRCoeff > 0)
    {
        for (n = nLRCoeff - 1, p = pLRCoeff[nLRCoeff-1], L = 0.0; n >= 0; n--)
        {
            while (p > pLRCoeff[n]) 
            {
                L *= edR;
                dLdR *= edR;
                p--;
            }
            L -= LRCoeff[n];
            dLdR += LRCoeff[n] * double(p);
        }
        L *= pow(edR, pLRCoeff[0]);
        dLdR *= pow(edR, pLRCoeff[0] + 1);
        if (A_ex != 0)
        {
            z = -A_ex * pow(Ra, gamma) * exp(-beta * Ra);
            L += z;
            dLdR += (gamma * edR - beta) * z;
        }
    }
    L += Uinf;
    /*if (pnC > pLRCoeff[nLRCoeff - 1])
    {
        pLRCoeff[nLRCoeff] = pnC;
        LRCoeff[nLRCoeff++] = 0.0;
    }
    else 
    {*/
    for (n = nLRCoeff -1; (n>=0 ? !LRCFree[n] : false); n--) ;
    pnC = pLRCoeff[(n>=0 ? n : n = nLRCoeff - 1)];
    //}
    LRCoeff[n] += (D = (dAdR - dLdR) / (double(pnC) * (x = pow(Ra, -pnC - 1))));
    /*if (ADLRC)
    {
        double dADdR = 1e99, dRaq = 1.0 / (Ra * Ra), cM, ldADdR, dM = 0.1;
        do
        {
            ldADdR = dADdR;
            for (n = nLRCoeff; n <= nLRCoeff + 2; n++) LRCoeff[n] = pow(LRCoeff[n-1] / LRCoeff[n-2], 3) * LRCoeff[n-3];
            for (n = nLRCoeff + 2, cM = 16, dADdR = 0.0; n >= nLRCoeff; n--, cM -= 2)
            {
                dADdR += cM * LRCoeff[n];
                dADdR *= dRaq;
            }
            if (dADdR > 5.0 * LRCoeff[nLRCoeff - 1])
            {
                printf("Error sdConnectLR: iterative procedure failed, no ADLRC!\n");
                ADLRC = false;
                break;
            }
            LRCoeff[nLRCoeff-1] -= dADdR * dM;
            L += x * dADdR * Ra * dM;
            //printf("dADdR=%g\n", dADdR);
        } 
        while (ldADdR - dADdR  > 1.0);
        if (ADLRC) 
        {
            for (n = nLRCoeff; n <= nLRCoeff + 2; n++) pLRCoeff[n] = pLRCoeff[n-1] + 2;
            for (n = nLRCoeff + 2, dADdR = 0.0; n >= nLRCoeff; n--)
            {
                dADdR += LRCoeff[n];
                dADdR *= dRaq;
                LRCFree[n] = false;
            }
            Tm -= x * dADdR;
            nLRCoeff += 3;
        }
    }*/
    Tm += (L - A - D * Ra * x);
    Lock->unlock();
}

void AnaPot::cdConnectLR(int pnC)
{
    int n, p, m=0;
    double dxdR, z, A = 0.0, dAdR = 0.0, x, L = 0.0, dLdR = 0.0, edR, D, F;
    Lock->lock();
    z = Ra + b * Rm;
    x = (Ra - Rm) / z;
    dxdR = (1.0 - x) / z;
    for (n = nPotCoeff - 1; n >= 0; n--)
    {
        dAdR += double(n+1) * PotCoeff[n];
        if (n>0) dAdR *= x;
        A += PotCoeff[n];
        A *= x;
    }
    dAdR *= dxdR;
    A += Tm;
    edR = 1.0 / Ra;
    if (nLRCoeff > 0)
    {
        if (pnC < pLRCoeff[0]) pnC = pLRCoeff[nLRCoeff - 1] + 2;
        else for (m = nLRCoeff - 1; pnC < pLRCoeff[m]; m--) ;
        for (n=m-1, p = pLRCoeff[n], L = 0.0; n >= 0; n--)
        {
            while (p > pLRCoeff[n]) 
            {
                L *= edR;
                dLdR *= edR;
                p--;
            }
            L -= LRCoeff[n];
            dLdR += LRCoeff[n] * double(p);
        }
        L *= pow(edR, pLRCoeff[0]);
        dLdR *= pow(edR, pLRCoeff[0] + 1);
    }
    else if (pnC < 6) pnC = 6;
    if (m > nLRCoeff - 2)
    {
        if (nLRCoeff > 0)
        {
            int *PB = pLRCoeff;
            double *CB = LRCoeff;
            pLRCoeff = new int[m+2];
            LRCoeff = new double[m+2];
            for (n=0; n < nLRCoeff; n++)
            {
                pLRCoeff[n] = PB[n];
                LRCoeff[n] = CB[n];
            }
            delete[] PB;
            delete[] CB;
        }
        else
        {
            pLRCoeff = new int[2];
            LRCoeff = new double[2];
        }
    }
    L += Uinf;
    pLRCoeff[m] = pnC;
    pLRCoeff[m + 1] = pnC + 2;
    x = pow(Ra, pnC + 1);
    D = dAdR - dLdR;
    F = D + double(pnC) / Ra * (A-L);
    //printf("Ra=%e, pnC=%d, x=%g, D=%g, F=%g\n", Ra, pnC, x, D, F);
    LRCoeff[m] = x / double(pnC) * (D - 0.5 * double(pnC + 2) * F);
    LRCoeff[m + 1] = 0.5 * Ra * Ra * x * F;
    nLRCoeff = m + 2;
    Lock->unlock();
}

void AnaPot::cdConnectSR(const bool)
{
    if (nPotCoeff == 0) return;
    double z, iR, d;
    double X, x;
    double dXdR;
    int n;
    QMutexLocker Locker(Lock);
    for (d = 0.0; d >= 0.0; Ri += 0.01)
    {
        if (Ri >= Ra) return;
        z = Ri + b * Rm;
        X = 1.0;
        x = (Ri - Rm) / z;
        dXdR = (1.0 - x) / z;
        for (n=0, d = 0.0; n < nPotCoeff; n++, X*=x) d += double(n+1) * PotCoeff[n] * X * dXdR;
    }
    Ri -= 0.01;
    //printf("Ri=%g, d=%g\n", Ri, d);
    iR = pow(Ri, -1.0 - iExp);
    iCoeff = d / (-iExp * iR);
    iOff = Point(Ri, 2) - iCoeff * iR * Ri;
}

void AnaPot::sConnect()
{
    Lock->lock();
    Tm += Point(Ra, 3) - Point(Ra, 2);
    iOff += Point(Ri, 2) - Point(Ri, 1);
    Lock->unlock();
}

void AnaPot::shift(double E)
{
    Uinf += E;
    iOff += E;
    Tm += E;
}

double AnaPot::getUinf()
{
    return Uinf;
}

int AnaPot::cAdRows()
{
    int n, AdRows;
    for (n = AdRows = 0; n < nLRCoeff; n++) if (LRCFree[n] && pLRCoeff[n] >= 12) AdRows++;
    return AdRows;
}

void AnaPot::calc_hCi_cP_EQS(int NEQ, double **EQS, double *Sig, double *EDiff, double hCi_cP)
{
    int n, m, k, l, i, j, N = nPotCoeff + nLRCoeff;
    double CC, F, z, x, X, dxdR;
    z = Ra + b * Rm;
    x = (Ra - Rm) / z;
    dxdR = (1.0 - x) / z;
    for (n=i=0; i < nLRCoeff; i++) if (pLRCoeff[i] >= 12 && LRCFree[i]) n++; 
    for (m=i=0; i < nLRCoeff; i++) if (pLRCoeff[i] >= 12 && LRCFree[i])
    {
        Sig[j = NEQ - n + m++] = 1.0 / ((CC = pow(LRCoeff[i-1] / LRCoeff[i-2], 3.0) * LRCoeff[i-3]) * hCi_cP);
        EDiff[j] = CC - LRCoeff[i];
        if (m==n)
        {
            F = pow(Ra, pLRCoeff[i] + 1) / double(pLRCoeff[i]);
            for (k=l=0, X = F * dxdR; k < nPotCoeff; k++, X*=x) if (SMap[l] == k) EQS[j][l++] = double(k+1) * X;
            for (k=0, X = F * (x = 1.0 / Ra); (l < NumFreePar ? SMap[l] < N : false); k++, X*=x) 
                if (pLRCoeff[SMap[l] - nPotCoeff] == k) EQS[j][l++] = double(-k) * X;
            if (ExcIntFree) EQS[j][l++] = (gamma * x - beta) * pow(Ra, gamma) * exp(-beta * Ra);
        }
        else l=0;
        while (l < NumFreePar) EQS[j][l++] = 0.0;
        for (l=0; (l < NumFreePar ? SMap[l] - nPotCoeff < i-3 : false); l++) ;
        if (SMap[l] - nPotCoeff == i-3) EQS[j][l++] -= CC / LRCoeff[i-3];
        if (SMap[l] - nPotCoeff == i-2) EQS[j][l++] += 3.0 * CC / LRCoeff[i-2];
        if (SMap[l] - nPotCoeff == i-1) EQS[j][l++] -= 3.0 * CC / LRCoeff[i-1];
        if (SMap[l] - nPotCoeff == i) EQS[j][l] = 1.0;
    }
}

double AnaPot::calcPotCoeff(int N, double *R, double *U, double *Sig, int n, double nRm, double nb)
{
    int i;
    double *UC = new double[N], FQSv, FQSn, TmC, PotCoeffC[n];
    b = nb;
    Rm = nRm;
    Lock->lock();
    if (Ri == 0.0) Ri = R[0];
    if (Ra == 0.0) Ra = R[N-1];
    if (nPotCoeff > 0) delete[] PotCoeff;
    PotCoeff = new double[nPotCoeff = n];
    FQSv = 2.0 * (FQSn = fitPotCoeff(N, R, U, n, Rm, b, PotCoeff, Tm, Sig));
    while (FQSn < 0.9999 * FQSv)
    {
        FQSv = FQSn;
        for (i=0; i<N; i++) UC[i] = U[i] - Point(R[i], 2);
        FQSn = fitPotCoeff(N, R, UC, n, Rm, b, PotCoeffC, TmC, Sig);
        for (i=0, Tm += TmC; i<n; i++) PotCoeff[i] += PotCoeffC[i];
    }
    delete[] UC;
    Lock->unlock();
    return FQSn;
    
    /*int i, j, k, m=2*n+1;
    double x, F[m], **M;
    M = Create(n + 1, n + 2);

    for (i=0; i<=n; i++) for (j=0; j<n+2; j++) M[i][j] = 0.0;
    for (i=0, F[0]=1.0; i<N; i++)
    {
        x = (R[i] - Rm) / (R[i] + b * Rm);
        for (j=1; j<m; j++) F[j] = F[j-1] * x;
        for (j=0; j<=n; j++) 
        {
            for (k=0; k<=n; k++) M[j][k] += F[k+j];
            M[j][n+1] += F[j] * U[i];
        }
    }
    SolvLinEqS(M, n+1, 1e-15);
    Tm = M[0][n+1];
    for (i=0; i<n; i++) PotCoeff[i] = M[i+1][n+1];
    Destroy(M, n + 1);*/
    
    //printf("Tm=%f\n", Tm);
    //for (i=0; i<n; i++) printf("a[%d]=%f\n", i, PotCoeff[i]);
}

void AnaPot::fitAdCorrToRmb()
{
    if (NumFreeAdCorr > 0)
    {
        int NPoints = 10000, n, m, a, SStart = nPotCoeff + nLRCoeff - (ExcIntFree ? 0 : 1), SMap[NAdCorr];
        double RMax = 20.0, *adCorrC = new double[NPoints], **EQS = Create(NPoints, NAdCorr), C[NAdCorr];
        double *Sig = new double[NPoints], *AD = new double[NPoints], FQS, oFQS;
        for (n=a=0, m = SStart; n < NAdCorr; n++, m++) if (adCorrFree[n]) SMap[a++] = m;
        for (n=0; n < NPoints; n++) 
        {
            adCorrC[n] = 0.0;
            Sig[n] = 0.01;
        }
        calcAdCorrPot(adCorrC, adCorrC, Ri, (RMax - Ri) / double(NPoints - 1), NPoints, 1.0, 1.0, NAdCorr, 1, PAdCorr, 0.0, 0.0, 
                      AdCorr_Rm, AdCorr_b, adCorr);
        if (AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
        calcS(NPoints, Ri, RMax);
        for (a=0; a < NAdCorr; a++) C[a] = adCorr[a] / AdCorrWF[a];
        for (n=0; n < NPoints; n++) 
        {
            for (a=0; a < NAdCorr; a++) if (!adCorrFree[a]) adCorrC[n] -= C[a] * S[SStart + a][n];
            for (a=0; a < NumFreeAdCorr; a++) EQS[n][a] = S[SMap[a]][n];
        }
        Destroy(S, numFreeCoeff);
        S=0;
        SvdFit(EQS, C, adCorrC, Sig, NumFreeAdCorr, NPoints, 1e-12);
        for (a=0; a < NumFreeAdCorr; a++) adCorr[SMap[a] - SStart] = C[a];
        for (n=0, FQS = 0.0; n < NPoints; n++)
        {
            for (a=0, AD[n] = adCorrC[n]; a < NumFreeAdCorr; a++) AD[n] -= C[a] * EQS[n][a];
            FQS += AD[n] * AD[n];
        }
        printf("AnaPot::fitAdCorrToRmb(): FQS=%f\n", FQS *= 1e4);
        oFQS = 2.0 * FQS;
        while (FQS < 0.9999 *oFQS && FQS > 1.0)
        {
            oFQS = FQS;
            SvdFit(EQS, C, AD, Sig, NAdCorr, NPoints, 1e-12);
            for (a=0; a < NumFreeAdCorr; a++) adCorr[SMap[a] - SStart] += C[a];
            for (n=0, FQS = 0.0; n < NPoints; n++)
            {
                for (a=0, AD[n] = adCorrC[n]; a < NumFreeAdCorr; a++) AD[n] -= adCorr[SMap[a] - SStart] * EQS[n][a];
                FQS += AD[n] * AD[n];
            }
            printf("AnaPot::fitAdCorrToRmb(): FQS=%f\n", FQS *= 1e4);
        }
        for (a=0; a < NumFreeAdCorr; a++) adCorr[SMap[a] - SStart] *= AdCorrWF[SMap[a] - SStart];
        delete[] adCorrC;
        delete[] Sig;
        delete[] AD;
        Destroy(EQS, NPoints);
    }
    AdCorr_b = b;
    AdCorr_Rm = Rm;
}

double AnaPot::fitPotential(int NPoints, double R0, double h, double* U, double* Sig, int NCoeff, int NLRC, int* pLRC, 
                            double nRi, double nRa, double nRm, double nb, double nUinf, double niExp)
{
    int i, j, NC = NCoeff + NLRC - 1;
    double **EQS = Create(NPoints, NC), *UD = new double[NPoints];
    double oFQS = 1e99, FQS, C[NC], R;
    bool NewLRC = false;
    if (SMap != 0)
    {
        delete[] SMap;
        SMap = 0;
    }
    if (NCoeff != nPotCoeff)
    {
        if (PotCoeff != 0)
        {
            delete[] PotCoeff;
            delete[] CoeffFree;
        }
        PotCoeff = new double[nPotCoeff = NCoeff];
        CoeffFree = new bool[NCoeff];
    }
    if (nLRCoeff != NLRC)
    {
        if (LRCoeff != 0) 
        {
            delete[] LRCoeff;
            delete[] pLRCoeff;
            delete[] LRCFree;
        }
        LRCoeff = new double[nLRCoeff = NLRC];
        LRCFree = new bool[NLRC];
        pLRCoeff = pLRC;
        NewLRC = true;
    }
    else for (i=0; i < NLRC; i++) if (pLRCoeff[i] != pLRC[i]) 
    {
        NewLRC = true;
        pLRCoeff[i] = pLRC[i];
    }
    for (i=0; i < NCoeff; i++) 
    {
        PotCoeff[i] = 0.0;
        CoeffFree[i] = true;
    }
    if (NewLRC) 
    {
        for (i=0; i < NLRC; i++) 
        {
            LRCoeff[i] = 0.0;
            LRCFree[i] = true;
        }
        numFreeLRC = NLRC - 1;
        NumFreePar = numFreeCoeff = NCoeff + numFreeLRC;
        SMap = new int[numFreeCoeff];
        for (i=0; i < numFreeCoeff; i++) SMap[i] = i;
    }
    else 
    {
        SMap = new int[NCoeff + NLRC - 1];
        for (j=0; j < NCoeff; j++) SMap[j] = j;
        for (i = numFreeLRC = 0; i < NLRC - 1; i++) if (LRCFree[i]) 
        {
            numFreeLRC++;
            SMap[j++] = NCoeff + i;
            LRCoeff[i] = 0.0;
        }
        NumFreePar = numFreeCoeff = NCoeff + numFreeLRC;
    }
    Ri = nRi;
    Ra = nRa;
    Rm = nRm;
    b = nb;
    Uinf = nUinf;
    iExp = niExp;
    calcS(NPoints, R0, R0 + h * double(NPoints - 1));
    double WLRC[NLRC - 1];
    for (i=0; i < NLRC - 1; i++) if (!LRCFree[i]) WLRC[i] = LRCoeff[i] / SWF[NCoeff + i];
    for (i=0; i < NPoints; i++) 
    {
        UD[i] = U[i];
        for (j=0; j < numFreeCoeff; j++) EQS[i][j] = S[SMap[j]][i];
        for (j=0; j < NLRC - 1; j++) if (!LRCFree[j]) UD[i] -= S[NCoeff + j][i] * WLRC[j];
    }
    Destroy(S, NC);
    S=0;
    SvdFit(EQS, C, UD, Sig, numFreeCoeff, NPoints, 1e-12);
    updatePotential(C);
    for (i=0, FQS = 0.0, R = R0; i < NPoints; i++, R+=h)
    {
        UD[i] = U[i] - Point(R);
        FQS += UD[i] * UD[i] / (Sig[i] * Sig[i]);
    }
    printf("FQS = %f\n", FQS);
    //printf("C6=%f\n", C[20] * SWF[20]);
    while (FQS < 0.9999 * oFQS)
    {
        SvdFit(EQS, C, UD, Sig, numFreeCoeff, NPoints, 1e-12);
        oFQS = FQS;
        updatePotential(C);
        for (i=0, FQS = 0.0, R = R0; i < NPoints; i++, R+=h)
        {
            UD[i] = U[i] - Point(R);
            FQS += UD[i] * UD[i] / (Sig[i] * Sig[i]);
        }
        printf("FQS = %f\n", FQS);
    }
    //printf("C6=%f\n", C[20] * SWF[20]);
    delete[] UD;
    Destroy(EQS, NPoints);
    printf("b=%f, FQS=%f\n", b, FQS);
    return FQS;
}

void AnaPot::calcS(int numPoints, double Rmin, double Rmax)
{
    int n, m, k, NOC = nPotCoeff + nLRCoeff - 1; 
    if (ExcIntFree) NOC++;
    NumS = NFitCoeff = NOC + NAdCorr;
    for (n = nLRCoeff - 1; (n>=0 ? !LRCFree[n] : false); n--) ;
    int PLLRC = pLRCoeff[(n>=0 ? n : nLRCoeff - 1)];
    double R, AO[nPotCoeff], AF[nPotCoeff], IF[nPotCoeff], IO[nPotCoeff], LRCO[nLRCoeff - 1], LRCF[nLRCoeff - 1];
    double LLRCO = -pow(Ra, -PLLRC), x = 0.0, xF = 0.0, DRa = 1.0 / (Ra + b * Rm), XA = (Ra - Rm) * DRa, XFA = 1.0;
    double dRi = 1.0 / (Ri + b * Rm), XI = (Ri - Rm) * dRi, XFI = 1.0, LLRCd = -Ra / (double(PLLRC) * LLRCO);
    double ICO = pow(Ri, -iExp), dICd = -Ri / (iExp * ICO), dRa = 1.0 / Ra, dXIdRi = (1.0 - XI) * dRi;
    double dXAdRa = (1.0 - XA) * DRa, ExcIntF = 0.0, ExcIntO = 0.0, h = (Rmax - Rmin) / double(numPoints - 1);
    if (SWF != 0) delete[] SWF;
    SWF = new double[NOC];
    S = Create(NFitCoeff, numPoints);
    for (n=0; n < nPotCoeff; n++)
    {
        AF[n] = double(n+1) * XFA * dXAdRa * LLRCd;
        AO[n] = AF[n] * LLRCO - (XFA *= XA);
        IF[n] = XFI * dICd * double(n+1) * dXIdRi;
        IO[n] = AO[n] + (XFI *= XI) - ICO * IF[n];
    }
    for (n=k=0, R=1.0; n < nLRCoeff - 1; n++, k++)
    {
        while (k < pLRCoeff[n])
        {
            R *= dRa;
            k++;
        }
        LRCO[n] = -R;
        LRCF[n] = -double(k) * (R *= dRa) * LLRCd;
        LRCO[n] += (LRCF[n] * LLRCO);
    }
    if (ExcIntFree)
    {
        ExcIntO = -pow(Ra, gamma) * exp(-beta * Ra);
        ExcIntF = -(gamma * dRa - beta) * ExcIntO * LLRCd;
        ExcIntO += ExcIntF * LLRCO;
    }
    for (n=0; n < NOC; n++) SWF[n] = 0.0;
    for (n=0; n < NAdCorr; n++) AdCorrWF[n] = 0.0;
    for (n=0, R = Rmin; n < numPoints; n++, R+=h)
    {
        if (R <= Ri)
        {
            x = pow(R, -iExp);
            for (k=0; k < nPotCoeff; k++) SWF[k] += fabs(S[k][n] = x * IF[k] + IO[k]);
            for (m=0; m < nLRCoeff - 1; m++, k++) SWF[k] += fabs(S[k][n] = LRCO[m]);
        }
        else if (R < Ra)
        {
            x = (R - Rm) / (R + b * Rm);
            for (k=0, xF = x; k < nPotCoeff; k++, xF *= x) SWF[k] += fabs(S[k][n] = AO[k] + xF);
            for (m=0; m < nLRCoeff - 1; m++, k++) SWF[k] += fabs(S[k][n] = LRCO[m]);
        }
        else
        {
            for (m=1, S[nPotCoeff][n] = -pow(DRa = 1.0 / R, k = pLRCoeff[0]); m < nLRCoeff - 1; m++)
                for (S[nPotCoeff + m][n] = S[nPotCoeff + m - 1][n]; k < pLRCoeff[m]; k++) 
                    S[nPotCoeff + m][n] *= DRa;
            for (xF = S[nPotCoeff + m - 1][n]; k < pLRCoeff[m]; k++) xF *= DRa;
            for (k=0; k < nPotCoeff; k++) SWF[k] += fabs(S[k][n] = AF[k] * xF);
            for (m=0; m < nLRCoeff - 1; m++, k++) SWF[k] += fabs(S[k][n] += (LRCF[m] * xF));
        }
        if (ExcIntFree)
        {
            if (R < Ra) S[k][n] = ExcIntO;
            else S[k][n] = -pow(R, gamma) * exp(-beta * R) + ExcIntF * xF;
            SWF[k] += fabs(S[k][n]);
            k++;
        }
        if (NAdCorr > 0)
        {
            if (NAdCorr > 1) x = (R - Rm) / (R + b * Rm);
            AdCorrWF[0] += fabs(S[k++][n] = pow(2.0 * Rm / (R + Rm), PAdCorr));
            for (m=1; m < NAdCorr; m++, k++) AdCorrWF[m] += fabs(S[k][n] = S[k-1][n] * x);
        }
    }
    for (n=0; n < NOC; n++) SWF[n] = 1.0 / SWF[n];
    for (n=0; n < NAdCorr; n++) AdCorrWF[n] = 1.0 / AdCorrWF[n];
    for (n=0; n < numPoints; n++)
    {
        for (k=0; k < NOC; k++) S[k][n] *= SWF[k];
        for (m=0; m < NAdCorr; m++, k++) S[k][n] *= AdCorrWF[m];
    }
}

void AnaPot::createSMap()
{
    int n, m;
    for (n = NumFreeAdCorr = 0; n < NAdCorr; n++) if (adCorrFree[n]) NumFreeAdCorr++;
    NumFreePar = numFreeCoeff + NumFreeAdCorr;
    NumFreeLRC = numFreeLRC;
    if (ExcIntFree)
    {
        NumFreePar++;
        NumFreeLRC++;
    }
    SMap = new int[NumFreePar];
    for (n=m=0; n < nPotCoeff; n++) if (CoeffFree[n]) SMap[m++] = n;
    for (n=0; n < nLRCoeff - 1; n++) if (LRCFree[n]) SMap[m++] = n + nPotCoeff;
    if (!LRCFree[n] && SMap[m-1] >= nPotCoeff) 
    {
        m--;
        NumFreeLRC--;
        NumFreePar--;
    }
    if (ExcIntFree) SMap[m++] = n + nPotCoeff;
    for (n=0; n < NAdCorr; n++) if (adCorrFree[n]) SMap[m++] = n + nPotCoeff + nLRCoeff - (ExcIntFree ? 0 : 1);
}

PotentialData* AnaPot::getPotentialData()
{
    int n;
    PotentialData *R = new PotentialData();
    R->PotType = analyticalPotential;
    Lock->lock();
    R->nPotCoeff = nPotCoeff;
    R->nLRCoeff = nLRCoeff;
    R->LRCoeff = new double[nLRCoeff];
    R->pLRCoeff = new int[nLRCoeff];
    R->LRCFree = new bool[nLRCoeff];
    for (n=0; n < nLRCoeff; n++)
    {
        R->LRCoeff[n] = LRCoeff[n];
        R->pLRCoeff[n] = pLRCoeff[n];
        R->LRCFree[n] = LRCFree[n];
    }
    R->NAdCorr = NAdCorr;
    R->TAdCorr = TAdCorr;
    R->PAdCorr = PAdCorr;
    R->NSpinRGamma = NSpinRGamma;
    R->Ri = Ri;
    R->Ra = Ra;
    R->iExp = iExp;
    R->iOff = iOff;
    R->iCoeff = iCoeff;
    R->b = b;
    R->Rm = Rm;
    R->Tm = Tm;
    R->De = Uinf;
    R->PotCoeff = new double[nPotCoeff];
    R->CFree = new bool[nPotCoeff];
    for (n=0; n < nPotCoeff; n++) 
    {
        R->PotCoeff[n] = PotCoeff[n];
        R->CFree[n] = CoeffFree[n];
    }
    R->A_ex = A_ex;
    R->beta = beta;
    R->gamma = gamma;
    R->q_e = q_e;
    R->q_f = q_f;
    if (NAdCorr > 0)
    {
        R->adCorr = new double[NAdCorr];
        R->adCorrFree = new bool[NAdCorr];
        for (n=0; n < NAdCorr; n++) 
        {
            R->adCorr[n] = adCorr[n];
            R->adCorrFree[n] = adCorrFree[n];
        }
    }
    R->RIso1AdCorr = RIso1AdCorr;
    R->RIso2AdCorr = RIso2AdCorr;
    R->SpinRGamma = (NSpinRGamma > 0 ? new double[NSpinRGamma] : 0);
    for (n=0; n < NSpinRGamma; n++) R->SpinRGamma[n] = SpinRGamma[n];
    R->SpinRR1 = SpinRR1;
    R->SpinRR2 = SpinRR2;
    R->EIFree = ExcIntFree;
    R->fitResult = Result;
    R->FQS = FQS;
    Lock->unlock();
    return R;
}

void AnaPot::restorePotential(double* bC)
{
    int n;
    Lock->lock();
    for (n=0; n < nPotCoeff; n++) PotCoeff[n] = bC[n];
    for (n=0; n < nLRCoeff - 1; n++) LRCoeff[n] = bC[nPotCoeff + n];
    for (n=0; n < NAdCorr; n++) adCorr[n] = bC[nPotCoeff + nLRCoeff - 1 + n];
    if (ExcIntFree) A_ex = bC[nPotCoeff + nLRCoeff + n];
    Lock->unlock();
    cdConnectLR1C();
    cdConnectSR(false);
}

void AnaPot::saveCoefficients(double*& bC)
{
    int n;
    if (bC == 0) bC = new double[nPotCoeff + nLRCoeff - (ExcIntFree ? 0 : 1) + NAdCorr];
    for (n=0; n < nPotCoeff; n++) bC[n] = PotCoeff[n];
    for (n=0; n < nLRCoeff - 1; n++) bC[nPotCoeff + n] = LRCoeff[n];
    for (n=0; n < NAdCorr; n++) bC[nPotCoeff + nLRCoeff + n - 1] = adCorr[n];
    if (ExcIntFree) bC[nPotCoeff + nLRCoeff + n] = A_ex;
}

void AnaPot::updatePotential(double *CD)
{
    Lock->lock();
    int n, m, N = NumFreePar, a;
    for (n=0; n<N; n++)
    {
        m = (SMap != 0 ? SMap[n] : n);
        if (m < nPotCoeff) PotCoeff[m] += SWF[m] * CD[n];
        else if (m < nPotCoeff + nLRCoeff - 1) LRCoeff[m - nPotCoeff] += SWF[m] * CD[n];
        else if (ExcIntFree && m < nPotCoeff + nLRCoeff) A_ex += SWF[m] * CD[n];
        else 
        {
            a = m - nPotCoeff - nLRCoeff + (ExcIntFree ? 0 : 1);
            adCorr[a] += AdCorrWF[a] * CD[n];
        }
    }
    Lock->unlock();
    cdConnectLR1C();
    cdConnectSR(false);
    //LRCoeff[nLRCoeff - 1] = 0.0;
}
