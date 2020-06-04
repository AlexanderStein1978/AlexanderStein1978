//
// C++ Implementation: TangToenniesPot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "TangToenniesPot.h"


TangToenniesPot::TangToenniesPot(PotFit *Fit) : PotWorker(Fit, TangToenniesPotential)
{
    Uinf = b = 0.0;
    C = SWF = 0;
    NC = 0;
    NA = 0;
    A = 0;
    A_ex = beta = gamma = 0.0;
    RMin = UMax = 0.0;
}

TangToenniesPot::TangToenniesPot(const TangToenniesPot& O) : PotWorker(O)
{
    int n;
    Uinf = O.Uinf;
    b = O.b;
    C = new double[NC = O.NC];
    for (n=0; n < NC; n++) C[n] = O.C[n];
    A = new double[NA = O.NA];
    for (n=0; n < NA; n++) A[n] = O.A[n];
    SWF = 0;
    A_ex = O.A_ex;
    beta = O.beta;
    gamma = O.gamma;
    RMin = O.RMin;
    UMax = O.UMax;
}

TangToenniesPot::~TangToenniesPot()
{
    if (C != 0) delete[] C;
    if (SWF != 0) delete[] SWF;
    if (A != 0) delete[] A;
}

void TangToenniesPot::AB(double B, double& Y, double& AA, double* D, double RM)
{
    //Translation of Fortran routine from K.T. Tang and J.P. Toennies, Z. Phys. D 1, 91 (1986)
    double F[NC], G[NC], SU, FR, GR, TK;
    int i, j, k;
    B *= RM;
    for (i=0; i < NC; i++)
    {
        SU = 0.0;
        FR = 1.0;
        for (j=1; j <= 2*(i+3); j++)
        {
            FR *= double(j);
            GR = pow(B, j) / FR;
            SU += GR;
        }
        G[i] = GR * exp(-B);
        F[i] = 1.0 - (1.0 + SU) * exp(-B);
    }
    AA = 0.0;
    Y = 0.0;
    for (k=0; k < NC; k++)
    {
        TK = 2*(k+3);
        AA += (F[k] * D[k] * TK / B - G[k] * D[k]);
        Y += (F[k] * D[k] * (TK / B - 1.0) - G[k] * D[k]);
    }
    AA *= exp(B);
    Y += 1.0;
}

bool TangToenniesPot::calcAb(double Re, double De)
{
    //Translation of Fortran routine from K.T. Tang and J.P. Toennies, Z. Phys. D 1, 91 (1986)
    int i, ik;
    double D[NC], DB, B1, Y1, Y2, A1, A2, R, YO, AO, BB, YY;
    for (i=0; i < NC; i++) D[i] = C[i] / (abs(De) * pow(Re, 2*(i+3)));
    b = 0.4;
    DB = 0.2;
    while(true)
    {
        for (ik = 1; ik <= 11; ik++)
        {
            B1 = b + DB;
            AB(b, Y1, A1, D, Re);
            AB(B1, Y2, A2, D, Re);
            R = Y2 / Y1;
            if (R <= 0.0) break;
            b += DB;
        }
        if (R > 0.0) return false;
        BB = b + DB / (1.0 - R);
        AB(BB, YO, AO, D, Re);
        YY = fabs(YO);
        if (YY <= 1e-12) break;
        DB /= 10.0;
    }
    A[0] = AO * abs(De);
    return true;
}


bool TangToenniesPot::calcPotCoeff(double* CC, int NCC, int FNC, double Re, double De)
{
    //De = 4.929e-3; Re = 8.828; CC[0] = 31.03e2; CC[1] = 37.92e4; CC[2] = 42.15e6;
    int i;
    if (NC != FNC)
    {
        if (C != 0) delete[] C;
        C = new double[NC = FNC];
    }
    if (A == 0) A = new double[1];
    for (i=0; i < NCC; i++) C[i] = CC[i];
    for ( ; i < NC; i++) C[i] = pow(C[i-1] / C[i-2], 3) * C[i-3];
    return calcAb(Re, De);
}

double TangToenniesPot::calcPotCoeff(int numPoints, double* R, double* U, double* Sig, int numCoeff, double nb)
{
    int n, m, NFC = numCoeff + (A_ex != 0.0 ? 2 : 1);
    b = nb;
    if (NC != numCoeff)
    {
        if (C != 0) delete[] C;
        C = new double[NC = numCoeff];
    }
    if (A == 0) A = new double[1];
    double **EQS = Create(numPoints, NFC), x[NFC], Calc, Diff, FQS;
    for (n=0; n <= numPoints; n++) getPreFakt(R[n], EQS[n]);
    SvdFit(EQS, x, U, Sig, NFC, numPoints, 1e-12);
    for (n=0, FQS = 0.0; n < numPoints; n++)
    {
        for (m=0, Calc = 0.0; m < NFC; m++) Calc += EQS[n][m] * x[m];
        Diff = U[n] - Calc;
        FQS += (Diff * Diff);
    }
    for (A[0] = x[n=0]; n < NC; n++) C[n] = x[n+1];
    if (A_ex != 0) A_ex = x[n+1];
    Destroy(EQS, numPoints);
    printf("FQS=%f\n", FQS);
    return FQS;
}

void TangToenniesPot::calcS(int numPoints, double Rmin, double Rmax)
{
    int n, R, r, m, NumCoefficients = getNFitCoefficients();
    double PreFakt[NumS = NumFreePar = NumCoefficients], x = 0.0, h = (Rmax - Rmin) / double(numPoints - 1);
    if (SWF != 0) delete[] SWF;
    SWF = new double[NumCoefficients];
    S = Create(NumCoefficients + NAdCorr, numPoints);
    for (n=0; n < NumCoefficients; n++) SWF[n] = 0.0;
    for (n=0; n < NAdCorr; n++) AdCorrWF[n] = 0.0;
    for (r=0, R = Rmin; r < numPoints; r++, R+=h)
    {
        getPreFakt(R, PreFakt);
        for (n=0; n < NumCoefficients; n++) SWF[n] += fabs(S[n][r] = PreFakt[n]);
        if (NAdCorr > 0)
        {
            if (NAdCorr > 1) x = (R - AdCorr_Rm) / (R + AdCorr_b * AdCorr_Rm);
            AdCorrWF[0] += fabs(S[n][r] = pow(2.0 * AdCorr_Rm / (R + AdCorr_Rm), PAdCorr));
            for (m=1; m < NAdCorr; m++, n++) AdCorrWF[m] += fabs(S[n][r] = S[n-1][r] * x);
        }
    }
    for (n=0; n < NumCoefficients; n++) SWF[n] = 1.0 / SWF[n];
    for (n=0; n < NAdCorr; n++) AdCorrWF[n] = 1.0 / AdCorrWF[n];
    for (r=0; r < numPoints; r++)
    {
        for (n=0; n < NumCoefficients; n++) S[n][r] *= SWF[n];
        for (m=0; m < NAdCorr; m++, n++) S[n][r] *= AdCorrWF[m];
    }
    NumCoefficients += NAdCorr;
}

void TangToenniesPot::getLRCoeff(double& R, int& numCoefficients, int*& Exponents, double*& Coefficients)
{
    int n;
    R = 0.0;
    Exponents = new int[numCoefficients = NC];
    Coefficients = C;
    for (n=1, Exponents[0] = 6; n < NC; n++) Exponents[n] = Exponents[n-1] + 2;
}

void TangToenniesPot::getMinimum(double& R, double& T)
{
    Lock->lock();
    double x1 = rmin, x2 = x1 + 1.0, x3 = x2 + 1.0;
    double y1 = Point(x1), y2 = Point(x2), y3 = Point(x3);
    int n;
    while (y3 < y2 && x3 < rmax)
    {
        x1 = x2;
        y1 = y2;
        x2 = x3;
        y2 = y3;
        y3 = Point(x3 = x2 + 1.0);
    }
    ParabInterpol(x1, y1, x2, y2, x3, y3, R, T);
    for (n=0; n < 1000 && fabs(x2 - R) >= 1e-8; n++)
    {
        if (R > x2)
        {
            x1 = x2;
            y1 = y2;
        }
        else 
        {
            x3 = x2;
            y3 = y2;
        }
        y2 = Point(x2 = R);
        ParabInterpol(x1, y1, x2, y2, x3, y3, R, T);
    }
    T = Point(R);
    Lock->unlock();
}

PotentialData* TangToenniesPot::getPotentialData()
{
    int n;
    PotentialData *R = new PotentialData();
    Lock->lock();
    R->A = new double[NA];
    for (n=0; n < NA; n++) R->A[n] = A[n];
    R->b = b;
    R->C = new double[NC];
    for (n=0; n < NC; n++) R->C[n] = C[n];
    R->Uinf = Uinf;
    if (NAdCorr > 0)
    {
        R->adCorr = new double[NAdCorr];
        for (n=0; n < NAdCorr; n++) R->adCorr[n] = adCorr[n];
    }
    R->RIso1AdCorr = RIso1AdCorr;
    R->RIso2AdCorr = RIso2AdCorr;
    R->AdCorr_b = AdCorr_b;
    R->AdCorr_Rm = AdCorr_Rm;
    R->N = NC;
    R->NA = NA;
    R->NAdCorr = NAdCorr;
    R->TAdCorr = TAdCorr;
    R->PAdCorr = PAdCorr;
    R->fitResult = Result;
    R->FQS = FQS;
    Lock->unlock();
    return R;
}

double TangToenniesPot::Point(double R, int)
{
    if (R <= RMin) return UMax;
    double U = 0.0, F[getNFitCoefficients()];
    int n, m;
    getPreFakt(R, F);
    for (n=0; n < NA; n++) U += A[n] * F[n];
    for (m=0; m < NC; m++) U += C[m] * F[n+m];
    if (A_ex != 0.0) U += A_ex * F[n+m];
    return U;
}

void TangToenniesPot::getPotCoeff(double *&rA, int &rNA, double& rb, double*& rC, int& N, double& rUinf)
{
    rA = A;
    rNA = NA;
    rb = b;
    rC = C;
    N = NC;
    rUinf = Uinf;
}

void TangToenniesPot::getPreFakt(double R, double* F)
{
    int n, m, i;
    double bR = b*R, Q = 1.0, dR = Q / R, dRq = dR * dR;
    for (i=1, F[0] = exp(-bR); i < NA; i++) F[i] = F[i-1] * dR; 
    for (n=0, m=1; n < NC; n++)
        for (F[n+i] = (n>0 ? F[n+i-1] : 1.0); m <= 2*(n+3); m++) F[n+i] += (Q *= (bR / m));
    for (n=i, Q = -pow(R, -6); n < NC + NA; n++, Q*= dRq) F[n] = Q * (1.0 - F[n] * F[0]);
    if (A_ex != 0) F[NC + NA] = pow(R, gamma) * exp(-beta * R);
}

double TangToenniesPot::improvePotential(int numPoints, double* R, double* U, double* Sig, double nb)
{
    int n, m, nc = getNFitCoefficients();
    double **EQS = Create(numPoints, nc), x[nc], Calc, aFQS = 1e99, FQS; 
    double *UD = new double[numPoints], Div, r1 = R[0], r2 = 0.0, r3, u1 = U[0], u2 = 0.0, u3, r, u;
    if (nb != 0.0) b = nb;
    for (n=0; n < NA; n++) x[n] = A[n];
    for (m=0; m < NC; m++) x[n+m] = C[m];
    if (A_ex != 0.0) x[n+m] = A_ex;
    for (n=0, FQS = 0.0; n < numPoints; n++)
    {
        getPreFakt(R[n], EQS[n]);
        for (m=0, UD[n] = U[n]; m < nc; m++) UD[n] -= EQS[n][m] * x[m];
        Div = UD[n] / Sig[n];
        FQS += Div * Div;
    }
    while (FQS < aFQS * 0.999)
    {
        printf("FQS=%g\n", aFQS = FQS);
        SvdFit(EQS, x, UD, Sig, nc, numPoints, 1e-12);
        for (n=0; n < NA; n++) x[n] = (A[n] += x[n]);
        for (m=0; m < NC; m++) x[n+m] = (C[m] += x[n+m]);
        if (A_ex != 0) x[n+m] = (A_ex += x[n+m]);
        for (n=0, FQS = 0.0; n < numPoints; n++)
        {
            for (m=0, Calc = 0.0; m < nc; m++) Calc += EQS[n][m] * x[m];
            Div = (UD[n] = (U[n] - Calc)) / Sig[n];
            FQS += (Div * Div);
        }
    }
    Destroy(EQS, numPoints);
    delete[] UD;
    printf("Finally: FQS=%g\n", FQS);
    do
    {
        r3 = r2;
        u3 = u2;
        r2 = r1;
        u2 = u1;
        r1 -= 0.1;
        u1 = Point(r1);
    }
    while (r1 >= rmin && u1 > u2);
    if (u1 < u2)
    {
        if (r2 == R[0]) u3 = Point(r3 = r2 + 0.1);
        do
        {
            ParabInterpol(r1, u1, r2, u2, r3, u3, r, u);
            if (r < r2)
            {
                r3 = r2;
                u3 = u2;
            }
            else
            {
                r1 = r2;
                u1 = u2;
            }
            r2 = r;
            u2 = Point(r);
        }
        while (r3 - r2 > 1e-12);
        RMin = r3;
        UMax = u3;
    }
    else RMin = UMax = 0.0;
    return FQS;
}

void TangToenniesPot::setPotCoeff(double *nA, int nNA, double nb, double* nC, int N, double nUinf)
{
    if (C != nC)
    {
        if (C != 0) delete[] C;
        C = nC;
    }
    if (A != nA)
    {
        if (A != 0) delete[] A;
        A = nA;
    }
    NA = nNA;
    b = nb;
    NC = N;
    Uinf = nUinf;
}

void TangToenniesPot::updateCoefficients(double* CD)
{
    int n, m; 
    for (n=0; n < NA; n++) A[n] += (SWF[n] * CD[n]);
    for (m=0; m < NC; m++) C[m] += (SWF[n+m] * CD[n+m]);
    if (A_ex != 0.0) A_ex += (SWF[n+m] + CD[n+m]);
}
