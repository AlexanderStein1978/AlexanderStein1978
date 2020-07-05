//
// C++ Implementation: PotWorker
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "potworker.h"
#include "SplinePoint.h"
#include "tableline.h"
#include "utils.h"
#include "isotab.h"
#include "AnaPot.h"
#include "splinepot.h"
#include "TangToenniesPot.h"
#include "MTTPot.h"
#include "MLRPot.h"
#include "PotFit.h"
#include "tools.h"
#include "potentialdata.h"

#include <QMutex>


PotWorker::PotWorker(PotFit *fit, PotentialType PotType)
{
    Type = PotType;
    SMap = 0;
    NumFreePar = 0;
    NumFreeAdCorr = NumFreeLRC = 0;
    Fit = fit;
    NSpinRGamma = 0;
    NAdCorr = PAdCorr = TAdCorr = 0;
    NFitCoeff = 0;
    AdCorrWF = 0;
    SpinRb = SpinRR1 = SpinRR2 = SpinRRm = 0.0;
    SpinRGamma = 0;
    RIso1AdCorr = RIso2AdCorr = AdCorr_b = AdCorr_Rm = 0.0;
    adCorr = 0;
    adCorrFree = 0;
    S = 0;
    Maxv = MaxJ = 0;
    nEs = Mv = nGS = 0;
    mJ = NEL = nGP = nGL = 0;
     PS = mv = 0;
    EL = 0;
    GL = 0;
    q_e = q_f = 0.0;
    Omega = 0;
    NIso = 0;
    el_S = 0;
    SFQSU = 0;
    BadListData = 0;
    BadListNR = 0;
    ShowBadList = false;
    points = 0;
    L = 0;
    NLRC = 0;
    PLRC = 0;
    LRC = 0;
    Uinf = 0.0;
    iO = iA = 0.0;
    IsoT = 0;
    numSplinePoints = 0;
    iExp = 6.0;
    Ro = 0.0;
    shortRVariable = true;
    LRCfree = 0;
    Lock = new QMutex;
    BLLock = new QMutex;
    Ndev = 0;
    Row = 0;
    dev = RWErr = 0;
    WFS = WWFS = 0;
    calcDiagFuncs = false;
    ES = 0;
    FQS = WeightSum = 0.0;
    NWeightF = 0;
    WeightFRow = 0;
    WeightF = 0;
    reCalcLC = reCalcSC = false;
    WantedValues = 0;
    NumWantedValues = 0;
}

PotWorker::PotWorker(const PotWorker& C) : QObject()
{
    int n;
    Type = C.Type; 
    SMap = 0;
    NumFreePar = 0;
    NumFreeAdCorr = NumFreeLRC = 0;
    NFitCoeff = 0;
    if ((NSpinRGamma = C.NSpinRGamma) > 0)
    {
        SpinRGamma = new double[NSpinRGamma];
        for (n=0; n < NSpinRGamma; n++) SpinRGamma[n] = C.SpinRGamma[n];
    }
    else SpinRGamma = 0;
    if ((NAdCorr = C.NAdCorr) > 0)
    {
        adCorr = new double[NAdCorr];
        adCorrFree = new bool[NAdCorr];
        for (n=0; n < NAdCorr; n++)
        {
            adCorr[n] = C.adCorr[n];
            adCorrFree[n] = C.adCorrFree[n];
        }
    }
    else
    {
        adCorr = 0;
        adCorrFree = 0;
    }
    Fit = 0;
    PAdCorr = C.PAdCorr;
    TAdCorr = C.TAdCorr;
    AdCorrWF = 0;
    SpinRb = C.SpinRb;
    SpinRR1 = C.SpinRR1;
    SpinRR2 = C.SpinRR2;
    SpinRRm = C.SpinRRm;
    RIso1AdCorr = C.RIso1AdCorr;
    RIso2AdCorr = C.RIso2AdCorr;
    AdCorr_b = C.AdCorr_b;
    AdCorr_Rm = C.AdCorr_Rm;
    nEs = Mv = nGS = 0;
    mJ = NEL = nGP = nGL = 0;
    EL = 0;
    GL = 0;
    PS = mv = 0;
    S = 0;
    q_e = C.q_e;
    q_f = C.q_f;
    Omega = C.Omega;
    NIso = 0;
    SFQSU = C.SFQSU;
    el_S = C.el_S;
    BadListData = 0;
    BadListNR = 0;
    ShowBadList = false;
    if ((numSplinePoints = C.numSplinePoints) > 0)
    {
        points = new SplinePoint[numSplinePoints];
        for (n=0; n < numSplinePoints; n++) points[n] = C.points[n];
    }
    else points = 0;
    if ((NLRC = C.NLRC) > 0)
    {
        PLRC = new int[NLRC];
        LRC = new double[NLRC];
        LRCfree = new bool[NLRC];
        for (n=0; n < NLRC; n++)
        {
            PLRC[n] = C.PLRC[n];
            LRC[n] = C.LRC[n];
            LRCfree[n] = C.LRCfree[n];
        }
    }
    else
    {
        PLRC = 0;
        LRC = 0;
        LRCfree = 0;
    }
    //printf("Mitte Potential\n");
    Uinf = C.Uinf;
    L = 0;
    iO = C.iO;
    iA = C.iA;
    iExp = C.iExp;
    IsoT = 0;
    Ro = C.Ro;
    shortRVariable = C.shortRVariable;
    Lock = new QMutex;
    BLLock = new QMutex;
    Ndev = 0;
    Row = 0;
    dev = RWErr = 0;
    WFS = WWFS = 0;
    calcDiagFuncs = false;
    ES = 0;
    FQS = WeightSum = 0.0;
    NWeightF = 0;
    WeightFRow = 0;
    WeightF = 0;
    reCalcSC = C.reCalcSC;
    reCalcLC = C.reCalcLC;
    WantedValues = 0;
    NumWantedValues = 0;
}

PotWorker::~PotWorker()
{
    if (NSpinRGamma > 0) delete[] SpinRGamma;
    if (NAdCorr > 0) 
    {
        delete[] adCorr;
        delete[] adCorrFree;
    }
    if (AdCorrWF != 0) delete[] AdCorrWF;
    if (S!=0) Destroy(S, NumS);
    if (SMap != 0) delete[] SMap;
    if (nEs != 0 || nGS != 0) destroyFitData();
    if (NIso > 0)
    {
        delete[] MU;
        delete[] IsoF;
    }
    if (BadListNR > 0) Destroy(BadListData, BadListNR);
    if (NLRC > 0)
    {
        delete[] PLRC;
        delete[] LRC;
        delete[] LRCfree;
    }
    if (points != 0) delete[] points;
    //printf("nach del aPot\n"); 
    if (IsoT != 0) delete IsoT;
    if (L!=0) Destroy(L, numSplinePoints);
    delete Lock;
    delete BLLock;
    if (dev != 0)
    {
        delete[] dev;
        delete[] Row;
        if (RWErr != 0) delete[] RWErr;
    }
    if (WFS != 0) delete[] WFS;
    if (WWFS != 0) delete[] WWFS;
    if (ES != 0) delete[] ES;
    if (NWeightF > 0)
    {
        delete[] WeightFRow;
        delete[] WeightF;
    }
    if (NumWantedValues > 0) delete[] WantedValues;
}

PotWorker* PotWorker::getCopy()
{
    switch (Type)
    {
        case NoPotential:
            return new PotWorker(*this);
        case analyticalPotential:
            return new AnaPot(*dynamic_cast<AnaPot*>(this));
        case SplinePotential:
            return new SplinePot(*dynamic_cast<SplinePot*>(this));
        case TangToenniesPotential:
            return new TangToenniesPot(*dynamic_cast<TangToenniesPot*>(this));
        case ModifiedTangToenniesPotential:
            return new MTTPot(*dynamic_cast<MTTPot*>(this));
        case MorseLongRangePotential:
            return new MLRPot(*dynamic_cast<MLRPot*>(this));
            break;
    }
    return 0;
}

bool PotWorker::setOmega(int nOmega, double S)
{
    if (Fit->isRunning()) return false;
    Omega = nOmega;
    el_S = S;
    return true;
}

bool PotWorker::setUinf(double nUinf)
{
    if (Fit->isRunning()) return false;
    Uinf = nUinf;
    return true;
}

void PotWorker::stopFit(bool Stop)
{
    if (Fit->isRunning() || StopFit) StopFit = Stop;
}

void PotWorker::getSplinePotForReading(int& rNSplinePoints, SplinePoint*& rpoints, int& rNLRC, int*& rpLRC, double*& rLRC, double& riA, 
                                       double& riO, double& riExp)
{
    printf("Lock->lock() PotWorker::getSplinePotForReading\n");
    Lock->lock();
    rNSplinePoints = numSplinePoints;
    rpoints = points;
    rNLRC = NLRC;
    rpLRC = PLRC;
    rLRC = LRC;
    riA = iA;
    riO = iO;
    riExp = iExp;
}

void PotWorker::getSplinePotForWriting(int& rNSplinePoints, SplinePoint*& rpoints, int& rNLRC, int*& rpLRC, double*& rLRC, double& riA,
                                       double& riO, double& riExp)
{
    if (Fit->isRunning())
    {
        rNSplinePoints = rNLRC = 0;
        rpLRC = 0;
        rLRC = 0;
        rpoints = 0;
        return;
    }
    rNSplinePoints = numSplinePoints;
    rpoints = points;
    rNLRC = NLRC;
    rpLRC = PLRC;
    rLRC = LRC;
    riA = iA;
    riO = iO;
    riExp = iExp;
}

bool PotWorker::setSplinePot(int nNSplinePoints, SplinePoint* npoints, int nNLRC, int* npLRC, double* nLRC, double niA, double niO, 
                             double niExp, double LongRangeF)
{
    if (Fit->isRunning()) return false;
    int n, i;
    if (points != npoints && npoints != 0) DestroyPotentials();
    numSplinePoints = nNSplinePoints;
    points = npoints;
    NLRC = nNLRC;
    if (nLRC != LRC)
    {
        if (PLRC != 0)
        {
            delete[] PLRC;
            delete[] LRC;
            delete[] LRCfree;
        }
        LRCfree = new bool[NLRC];
        for (n=0; n < NLRC; n++) LRCfree[n] = true;
    }
    PLRC = npLRC;
    LRC = nLRC;
    iA = niA;
    iO = niO;
    iExp = niExp;
    if (LongRangeF != 0.0)
    {
        calcLMatrix();
        for (n=0; n < numSplinePoints; n++)
        {
            for (i=0, points[n].yss = iA * L[n][0]; i < numSplinePoints ; i++) points[n].yss += points[i].y * L[n][i+1];
            points[n].yss += LongRangeF * L[n][numSplinePoints + 1];
            points[n].variable = true;
        }
    }
    return true;
}

AnaPot* PotWorker::getAnaPotForReading()
{
    printf("Lock->lock() PotWorker::getAnaPotForReading\n");
    if (Type == analyticalPotential) Lock->lock();
    return dynamic_cast<AnaPot*>(this);
}

AnaPot* PotWorker::getAnaPotForWriting()
{
    if (Fit->isRunning()) return 0;
    return dynamic_cast<AnaPot*>(this);
}

void PotWorker::unlock()
{
    printf("PotWorker::unlock\n");
    Lock->unlock();
}

double* PotWorker::getMCSE(int& N, int*& Rows)
{
    int s, n;
    for (s=N=0; s < nEs; s++) N += NEL[s];
    for (s=0; s < nGS; s++) N += nGL[s];
    double *R = new double[N];
    Rows = new int[N];
    for (s=N=0; s < nEs; s++) for (n=0; n < NEL[s]; n++)
    {
        Rows[N] = EL[s][n].Row;
        R[N++] = EL[s][n].E;
    }
    for (s=0; s < nGS; s++) for (n=0; n < nGL[s]; n++)
    {
        Rows[N] = GL[s][n].Row;
        R[N++] = GL[s][n].WN;
    }
    return R;
}

void PotWorker::DestroyPotentials()
{
    if (points != 0)
    {
        delete[] points;
        points = 0;
        numSplinePoints = 0;
        if (NLRC != 0)
        {
            delete[] LRC;
            delete[] LRCfree;
            LRC = 0;
            LRCfree = 0;
            NLRC = 0;
        }
        DestroyLS();
    }
}

void PotWorker::getAdCorrForWriting(int& N, double*& a, int& T, int& P, double& Iso1, double& Iso2, double &Rm, double &b)
{
    if (Fit->isRunning())
    {
        N = 0;
        return;
    }
    N = NAdCorr;
    a = adCorr;
    T = TAdCorr;
    P = PAdCorr;
    Iso1 = RIso1AdCorr;
    Iso2 = RIso2AdCorr;
    Rm = AdCorr_Rm;
    b = AdCorr_b;
}

void PotWorker::getAdCorrForReading(int& N, double*& a, int& T, int& P, double& Iso1, double& Iso2, double &Rm, double &b, 
                                    bool **rAdCorrFree, bool alreadyLocked)
{
    int n;
    if (!alreadyLocked)
    {
        Lock->lock();
        a = new double[NAdCorr];
        for (n=0; n < NAdCorr; n++) a[n] = adCorr[n];
    }
    else a = adCorr;
    N = NAdCorr;
    T = TAdCorr;
    P = PAdCorr;
    Iso1 = RIso1AdCorr;
    Iso2 = RIso2AdCorr;
    Rm = AdCorr_Rm;
    b = AdCorr_b;
    if (!alreadyLocked)
    {
        if (rAdCorrFree != 0)
        {
            *rAdCorrFree = new bool[N];
            for (n=0; n<N; n++) (*rAdCorrFree)[n] = adCorrFree[n];
        }
        Lock->unlock();
    }
    else if (rAdCorrFree != 0) *rAdCorrFree = adCorrFree; 
}

bool PotWorker::setAdCorr(int N, double* a, int T, int P, double RIso1, double RIso2, double Rm, double b, bool *AdCorrFree)
{
    if (Fit->isRunning()) return false;
    int n;
    if (adCorr != 0) 
    {
        delete[] adCorr;
        delete[] adCorrFree;
    }
    if (AdCorrWF != 0)
    {
        delete[] AdCorrWF;
        AdCorrWF = 0;
    }
    NAdCorr = N;
    adCorr = a;
    TAdCorr = T;
    PAdCorr = P;
    RIso1AdCorr = RIso1;
    RIso2AdCorr = RIso2;
    if (AdCorrFree == 0)
    {
        adCorrFree = new bool[N];
        for (n=0; n<N; n++) adCorrFree[n] = true;
    }
    else adCorrFree = AdCorrFree;
    for (n = NumFreeAdCorr = 0; n<N; n++) if (AdCorrFree[n]) NumFreeAdCorr++;
    if (Type == analyticalPotential)
    {
        double nTm, nDe, *C;
        int N;
        AnaPot *APot = dynamic_cast<AnaPot*>(this);
        APot->getPotCoeff(AdCorr_b, AdCorr_Rm, nTm, nDe, N, C);
        if (b != 1e-30 && (fabs(AdCorr_b - b) > 1e-13 || fabs(AdCorr_Rm - Rm) > 1e-13)) APot->fitAdCorrToRmb();
    }
    else
    {
        AdCorr_b = b;
        AdCorr_Rm = Rm;
    }
    return true;
}

void PotWorker::getLRCoeffForReading(int& N, int*& pLRC, double*& LRC, double *Ra, bool alreadyLocked)
{
    if (!alreadyLocked) Lock->lock();
    int n, *PLRC;
    double *RLRC, R;
    getLRCoeff(R, N, PLRC, RLRC);
    if (alreadyLocked)
    {
        pLRC = PLRC;
        LRC = RLRC;
    }
    else
    {
        pLRC = new int[N];
        LRC = new double[N];
        for (n=0; n<N; n++) 
        {
            pLRC[n] = PLRC[n];
            LRC[n] = RLRC[n];
        }
        Lock->unlock();
        if (Type == TangToenniesPotential) delete[] PLRC;
    }
    if (Ra != 0) *Ra = R; 
}

void PotWorker::getLRCoeffForWriting(int& N, int*& pLRC, double*& LRC, bool** LRCFree, double *rRa)
{
    if (Fit->isRunning())
    {
        N = 0;
        return;
    }
    double R;
    getLRCoeff(R, N, pLRC, LRC);
    if (LRCFree != 0) *LRCFree = getLRCFree();
    if (rRa != 0) *rRa = R;
}

PotWorker* PotWorker::scalePotential(double, double)
{
    printf("Error: Scaling is not implemented for Potentials of type %d\n", Type);
    return 0;
}

void PotWorker::fastWaveFunc(double *U, double E, double PM, double h, int NP, double *WF, double &Ri, 
                              double &Ra, double &ISq, double &F, int &i, int &c, int &a)
{
    double Sq, Y0, Y1, Y2, eins = 1.0, ED, M = 12.0 * PM;
    int n;
    getMinMaxInt(U, E, M, NP, i, a);
    if (a >= NP - 1)
    {
        for (a = NP - 1; (a>i? U[a-1] > U[a] || U[a] < E : false); a--) ;
        if (a==i) 
        {
            printf("E=%f is too high!\n", E);
            c=i;
            return;
        }
    }
    //printf("i=%d, a=%d\n", i, a);
    Ra = Ri + double(a) * h;
    Ri += double(i) * h;
    Y1 = 1e-5;
    Y0 = Y1 * exp(-sqrt(M*(U[a]-E)));
    c=a;
    WF[c] = Y0 / (eins - PM * (U[c] - E));
    Sq = WF[c] * WF[c];
    c--;
    WF[c] = Y1 / (eins - PM * (ED = U[c] - E));
    Sq += WF[c] * WF[c];
    for (c--; c>i; c--)
    {
        Y2 = M * ED * WF[c+1] + 2.0 * Y1 - Y0;
        WF[c] = Y2 / (eins - PM * (ED = U[c] - E));
        if (WF[c] <= WF[c+1]) break;
        Sq += WF[c] * WF[c];
        Y0 = Y1;
        Y1 = Y2;
    }
    F = eins / WF[c];
    Y1 = 1e-5;
    Y0 = Y1 * exp(-sqrt(M*(U[i]-E)));
    WF[n=i] = Y0 / (eins - PM * (U[i] - E));
    ISq = WF[n] * WF[n];
    n++;
    WF[n] = Y1 / (eins - PM * (ED = U[n] - E));
    ISq += WF[n] * WF[n];
    for (n++; n<=c; n++)
    {
        //WF[i-1] = Y[i-1] / (eins - PM * (ED = U[i-1] - E));
        //Y[i] = M * ED * WF[i-1] + 2.0 * Y[i-1] - Y[i-2];
        
        Y2 = M * ED * WF[n-1] + 2.0 * Y1 - Y0;
        WF[n] = Y2 / (eins - PM * (ED = U[n] - E));
        ISq += WF[n] * WF[n];
        //if (n-i==1241) printf("WF[%d]=%g, ISq=%g\n", n, WF[n], ISq);
        Y0 = Y1;
        Y1 = Y2;
        //if (WF[n-2] < WF[n-1] && WF[n-1] >= WF[n]) printf("UMax at R=%f\n", Ri + double(n-1-i) * h);
    }
    F *= WF[c];
    ISq += F*F * Sq;
    //printf("WF[%d]=%g, ISq=%g\n", a, WF[a], ISq);
    //printf("Ende getFastWaveFunc\n");
}

void PotWorker::setFitData(int* nmJ, int** nmv, int nMv, int nnEs, int* nNEL, TermEnergy** nEL, int nnGS, int* nnGP, int* nnGL, 
                           int** nPS, TableLine** NGL, bool *nES, int nMaxJ, int nMaxv)
{
    if (Fit->isRunning()) return;
    if (nEs > 0 || nGS > 0) destroyFitData();
    mJ = nmJ;
    mv = nmv;
    Mv = nMv;
    nEs = nnEs;
    NEL = nNEL;
    EL = nEL;
    nGS = nnGS;
    nGP = nnGP;
    nGL = nnGL;
    PS = nPS;
    GL = NGL;
    MaxJ = nMaxJ;
    Maxv = nMaxv;
    if (ES != 0) delete[] ES;
    ES = nES;
}

void PotWorker::CreateEK(double ****&E, double *****&K, int *mJ, int **mv, int N)
{
    //printf("Potential::CreateEK\n");
    if (NIso == 0) return;
    int I, J, v, F, NFC = int(2.0 * el_S) + 1;
    E = new double***[NFC];
    K = new double****[NFC];
    for (F=0; F < NFC; F++)
    {
        E[F] = new double **[NIso];
        K[F] = new double ***[NIso];
        for (I=0; I < NIso; I++) if (mJ[I] > -1)
        {
            E[F][I] = new double *[mJ[I] + 1];
            K[F][I] = new double **[mJ[I] + 1];
            for (J=0; J <= mJ[I]; J++) if (mv[I][J] > -1)
            {
                E[F][I][J] = new double[mv[I][J] + 1];
                K[F][I][J] = new double*[mv[I][J] + 1];
                for (v=0; v <= mv[I][J]; v++) K[F][I][J][v] = new double[N];
            }
        }
    }
    //printf("Ende CreateEK\n");
}

void PotWorker::calcE(double ****E, int *mJ, int **mv, int Mv, int NumWFPoints, double ****MinErr, bool ****SFQSU, double SFQSRad)
{
    int I, J, nv, FC;
    int NFC = int(2.0 * el_S) + 1;
    bool AdCorr = (NAdCorr > 0);
    double h = (rmax - rmin) / (NumWFPoints - 1), *RPot = new double[NumWFPoints], *Y = new double[NumWFPoints];
    double *U = getPoints(rmin, rmax, NumWFPoints), MF, *AdCorrPot = (AdCorr ? new double[NumWFPoints] : 0);
    double **WF = Create(Mv + 1, NumWFPoints);
    //printf("Nach Create WF\n");
    //aPot->GetMin(R, M);
    double *F = new double[Mv+1];
    double *SQ = new double[Mv+1];
    int *CP = new int[Mv+1], SFQSP = (SFQSU != 0 ? int((SFQSRad - rmin) / h) : 0);
    int *ri = new int[Mv+1];
    int *ra = new int[Mv+1], v;
    //printf("Mv=%d, mv[6][35]=%d\n", Mv, mv[6][35]);
    for (FC = 0; FC < NFC; FC++) for (I=0; I < NIso; I++) 
    {
        MF = IsoF[I] * h*h;
        if (AdCorr) CalcAdCorrPot(U, AdCorrPot, rmin, h, NumWFPoints, IsoT->mIso1[I], IsoT->mIso2[I]);
        else AdCorrPot = U;
        for (J=0; J <= mJ[I]; J++) if (mv[I][J] > -1)
        {
            //if (mv[I][J] >= 60) *DebugStream << "Iso=" << QString::number(I) << " J=" << QString::number(J) << '\n'; 
            
            CRotPot(J, I, AdCorrPot, RPot, rmin, h, NumWFPoints, 0, (NFC > 1 ? FC + 1 : 0));
            NumerovCooley(Y, MF, nv = mv[I][J] + 1, E[FC][I][J], WF, ri, ra, CP, F, SQ, NumWFPoints, RPot,
                          (MinErr != 0 ? MinErr[FC][I][J] : 0));
            if (nv < mv[I][J] + 1)
            {
                printf("For J=%d and Iso%d v=%d is not existing, highest v is %d\n",
                       J, I, mv[I][J], nv - 1);
                for (v = nv; v <= mv[I][J]; v++) E[FC][I][J][v] = 0.0;
            }
            if (SFQSU != 0) 
            {
                for (v=0; v < nv; v++) SFQSU[FC][I][J][v] = (CP[v] > SFQSP);
                while (v <= mv[I][J]) SFQSU[FC][I][J][v++] = false;
            }
        }
    }
    //for (I=0; I <= mv[6][35]; I++) printf("E[6][35][%d]=%f\n", I, E[6][35][I]);
    //printf("Vor delete\n");
    delete[] F;
    delete[] U;
    Destroy(WF, Mv + 1);
    delete[] SQ;
    delete[] CP;
    delete[] ri;
    delete[] ra;
     if (AdCorr) delete[] AdCorrPot;
    delete[] RPot;
    delete[] Y;
    //printf("Ende\n");
}

void PotWorker::CalcAdCorrPot(double* InPot, double* OutPot, double Rmin, double Step, int NPoints, double MIso1, double MIso2)
{
    Lock->lock();
    calcAdCorrPot(InPot, OutPot, Rmin, Step, NPoints, MIso1, MIso2, NAdCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, 
                  AdCorr_Rm, AdCorr_b, adCorr);
    Lock->unlock();
}

int PotWorker::getRefIso()
{
    if (IsoT != 0) return IsoT->refIso;
    return 0;
}

void PotWorker::CreateFitEQS(int &N, double **&EQS, double *&EO, double *&DE, double *&Sig, double *&Err,
                                   int **&QN, int NP, int nEs, int *NEL, TermEnergy **EL, int nGS, 
                                      int *nGP, int *nGL, int **PS, TableLine **GL, int AdRows)
{
    //printf("Potential::CreateSplineFitEQS\n");
    int i, n, m, l, d, NFC = int(2.0 * el_S) + 1;
    //double F;
    for (N=n=0; n < nEs; n++) N += NEL[n];
    //printf("1. N=%d\n", N);
    for (n=0; n < nGS; n++) if (nGL[n] > 0)
    {
        N += nGL[n];
        /*for (p=0; p < nGP[n] - 1; p++)
        {
            l = PS[n][p+1] - PS[n][p];
            N += l * (l+1) / 2 - l;
        }
        l = nGL[n] - PS[n][p];
        N += l * (l+1) / 2 - l;*/
        //printf("n=%d, N=%d\n", n, N);
    }
    //printf("Nach erster Schleife, N=%d\n", N);
    EQS = Create(N + AdRows + NumWantedValues, NP);
    EO = new double[N];
    DE = new double[N + AdRows + NumWantedValues];
    QN = CreateInt(N, 6);
    Sig = new double[N + AdRows + NumWantedValues];
    Err = new double[N];
    for (n=m=0; n < nEs; n++) for (l=0; l < NEL[n]; l++)
    {
        EO[m] = EL[n][l].E;
        Sig[m] = 1.0 / EL[n][l].err;
        Err[m] = EL[n][l].err;
        QN[m][0] = EL[n][l].Iso;
        QN[m][1] = EL[n][l].J;
        QN[m][2] = EL[n][l].v;
        QN[m][5] = (EL[n][l].FC < NFC ? EL[n][l].FC : 0);
        QN[m++][3] = -1;
    }
    //printf("Vor zweiter\n");
    for (n=0; n < nGS; n++) for (l=0; l < nGP[n]; l++)
    {
        d = (l+1 < nGP[n] ? PS[n][l+1] : nGL[n]) - PS[n][l];
        //F = double(d) / (d*(d-1)/2);
        for (i = PS[n][l]; i < PS[n][l] + d; i++) 
            //for (j=i+1; j < PS[n][l] + d; j++)
        {
            EO[m] = GL[n][i].WN;// - GL[n][j].WN;
            //Err[m] = sqrt(GL[n][i].err * GL[n][i].err + GL[n][j].err * GL[n][j].err);
            //Sig[m] = F / Err[m];
            Sig[m] = 1.0 / (Err[m] = GL[n][i].err);
            QN[m][0] = GL[n][i].Iso;
            QN[m][1] = GL[n][i].Jss;//j
            QN[m][2] = GL[n][i].vss;//j
            QN[m][3] = l;
            QN[m][4] = d;
            //QN[m][3] = GL[n][i].Jss;
            //QN[m][4] = GL[n][i].vss;
            QN[m++][5] = (GL[n][i].FC < NFC ? GL[n][i].FC : 0);
        }
    }
    //printf("Ende CreateSplineFitEQS, m=%d\n", m);
}

void PotWorker::getSFQSU(bool****& rSFQSU, int*& rmJ, int**& rmv)
{
    rSFQSU = SFQSU;
    rmJ = mJ;
    rmv = mv;
}

double PotWorker::calcFQS(double ****E, int nES, int *NEL, TermEnergy **EL, int nGS, int *nGP, int *nGL,
                          int **PS, TableLine **GL, int NumWFPoints, bool calcDiag, int MJ, int Mv, double *****K,
                          int *rv, double *RSig, double TEF, double hCi_cP, bool ****SFQSU)
{
    //printf("calcFQS\n");
    double D, FQS = 0.0, SS, OT, Sig, CE;
    int n, m, l, M, NK, k, N, r=0, s=0, d, i, j, NFC = int(2.0 * el_S) + 1, FC;
    double *KS, KM = 0.0, ***W1 = 0, ***W2 = 0, F;
    SFQS = 0.0;
    BLLock->lock();
    if (!calcDiagFuncs) calcDiag = false;
    else if (calcDiag) 
    {
        W1 = Create(NIso, Mv + 1, MJ + 1);
        W2 = Create(NIso, Mv + 1, MJ + 1);
        for (n=0; n < NIso; n++) for (m=0; m <= Mv; m++) for (l=0; l <= MJ; l++)
                    W1[n][m][l] = W2[n][m][l] = 0.0;
    }
    if (rv != 0) 
    {
        KS = new double[NK = *rv];
        for (k=0; k < NK; k++) KS[k] = 0.0;
    }
    else
    {
        KS = 0;
        NK = 0;
    }
    //printf("Vor erster Schleife\n");
    for (n=N=0; n < nES; n++) N += NEL[n];
    for (n=0; n < nGS; n++) N += nGL[n];
    if (RWErr != 0 && (RSig == 0 || Ndev != N)) delete[] RWErr;
    if (RSig != 0 && (Ndev != N || RWErr == 0)) RWErr = new double[N];
    else if (RSig == 0) RWErr = 0;
    for (n=0; n < nES; n++) for (l=0; l < NEL[n]; l++) 
    {
        FC = (EL[n][l].FC >= 0 && EL[n][l].FC < NFC ? EL[n][l].FC : 0);
        EL[n][l].DevR = D = (EL[n][l].dev = EL[n][l].E - (CE = E[FC][EL[n][l].Iso][EL[n][l].J][EL[n][l].v])) / EL[n][l].err;
        if (CE != 0.0) 
        {
            if (RSig != 0)
            {
                RWErr[r] = sqrt(EL[n][l].err * EL[n][l].err + 0.3 * EL[n][l].dev * EL[n][l].dev);
                RSig[s] = 1.0 / RWErr[r++];
                FQS += EL[n][l].dev * EL[n][l].dev * RSig[s] * RSig[s] * TEF * EL[n][l].WeightFactor;
                s++;
            }
            else FQS += D*D * TEF * EL[n][l].WeightFactor;
            if (SFQSU != 0 ? SFQSU[FC][EL[n][l].Iso][EL[n][l].J][EL[n][l].v] : false) SFQS += D*D * EL[n][l].WeightFactor;
            if (D > 1E3) printf("Error: D=%g, Iso=%d, J=%d, v=%d, E=%f\n", D, EL[n][l].Iso, EL[n][l].J, EL[n][l].v, EL[n][l].E);
        }
        else 
        {
            EL[n][l].DevR = 0.0;
            if (RSig != 0)
            {
                RWErr[r++] = EL[n][l].err;
                RSig[s++] = 1.0 / EL[n][l].err;
            }
        }
        if (calcDiag)
        {
            W1[EL[n][l].Iso][EL[n][l].v][EL[n][l].J] += 1.0 / EL[n][l].err;
            W2[EL[n][l].Iso][EL[n][l].v][EL[n][l].J] += D;
        }
        for (k=0; k < NK; k++) 
            KS[k] += fabs(D * K[(EL[n][l].FC >= 0 && EL[n][l].FC < NFC ? EL[n][l].FC : 0)][EL[n][l].Iso][EL[n][l].J][EL[n][l].v][k]);
        //printf("v=%d, J=%d, Iso=%d, D=%f\n", EL[n][l].v, EL[n][l].J, EL[n][l].Iso, D);
    }
    //printf("Vor 2.\n");
    for (n=0; n < nGS; n++) for (l=0; l < nGP[n]; l++)
    {
        SS = OT = 0.0;
        M = (l+1 < nGP[n] ? PS[n][l+1] : nGL[n]);
        for (m = PS[n][l]; m < M; m++)
        {
            if ((CE = E[(GL[n][m].FC >= 0 && GL[n][m].FC < NFC ? GL[n][m].FC : 0)][GL[n][m].Iso][GL[n][m].Jss][GL[n][m].vss]) == 0.0) continue;
            SS += Sig = 1.0 / (GL[n][m].err * GL[n][m].err);
            OT += Sig * (CE + GL[n][m].WN);
            //if (l==4) printf("m=%d, Iso=%d, Jss=%d, vss=%d, WN=%f, err=%f, E=%f\n", m, GL[n][m].Iso,
                //                GL[n][m].Jss, GL[n][m].vss, GL[n][m].WN, GL[n][m].err,
                     //            E[GL[n][m].Iso][GL[n][m].Jss][GL[n][m].vss]);
        }
        OT /= SS;
        for (m = PS[n][l]; m < M; m++)
        {
            FC = (GL[n][m].FC >= 0 && GL[n][m].FC < NFC ? GL[n][m].FC : 0);
            D = (GL[n][m].dev = GL[n][m].WN - OT + (CE = E[FC][GL[n][m].Iso][GL[n][m].Jss][GL[n][m].vss])) / GL[n][m].err;
            if (RSig != 0) RWErr[r] = sqrt(GL[n][m].err * GL[n][m].err + 0.03 * GL[n][m].dev * GL[n][m].dev);
            if (CE == 0.0) 
            {
                GL[n][m].DevR = 0.0;
                if (RSig != 0) RWErr[r++] = GL[n][m].err;
                continue;
            }
            GL[n][m].DevR = D;
            //if (fabs(D) > 10.0)
                //printf("v=%d, J=%d, Iso=%d, D=%f\n", GL[n][m].vss, GL[n][m].Jss, GL[n][m].Iso, D);
            FQS += (RSig != 0 ? GL[n][m].dev * GL[n][l].dev / (RWErr[r] * RWErr[r]) : D*D) * GL[n][m].WeightFact;
            if (SFQSU != 0 ? SFQSU[FC][GL[n][m].Iso][GL[n][m].Jss][GL[n][m].vss] : false) SFQS += D * D * GL[n][m].WeightFact; 
            r++;
            if (calcDiag)
            {
                W1[GL[n][m].Iso][GL[n][m].vss][GL[n][m].Jss] += 1.0 / GL[n][m].err;
                W2[GL[n][m].Iso][GL[n][m].vss][GL[n][m].Jss] += D;
            }
            for (k=0; k < NK; k++) 
                KS[k] += fabs(D * K[(GL[n][m].FC >= 0 && GL[n][m].FC < NFC ? GL[n][m].FC : 0)][GL[n][m].Iso][GL[n][m].Jss][GL[n][m].vss][k]);
        }
        if (RSig != 0)
        {
            d = M - PS[n][l];
            F = double(d) / (d*(d-1)/2);
            for (i = PS[n][l]; i<M; i++) for (j=i+1; j<M; j++) 
                RSig[s++] = F / sqrt(RWErr[r-M+i] * RWErr[r-M+i] + RWErr[r-M+j] * RWErr[r-M+j]);
        }    
        //printf("nGS=%d, l=%d, nGP[%d]=%d, FQS=%f, SS=%f, OT=%f\n", nGS, l, n, nGP[n], FQS, SS, OT);
        //printf("PS[n][l]=%d, M=%d\n", PS[n][l], M);
    }
    if (calcDiag)
    {
        if (WFS == 0) WFS = new double[NumWFPoints];
        if (WWFS == 0) WWFS = new double[NumWFPoints];
        double *WF = new double[NumWFPoints], *P = getPoints(rmin, rmax, NumWFPoints);
        double *U = new double[NumWFPoints], PM;
        for (n=0; n < NumWFPoints; n++) WFS[n] = WWFS[n] = 0.0;
        double h = (rmax - rmin) / (NumWFPoints - 1), M, Ri, Ra, ISq, F, F1, F2;
        int i, c, a;
        bool RC;
        for (n=0; n < NIso; n++)
        {
            M = IsoF[n] * h*h;
            PM = M * 0.0833333333333;
            for (m=0; m < MJ; m++)
            {
                RC = false;
                for (l=0; l < Mv; l++) if (W1[n][l][m] != 0.0)
                {
                    if (!RC) 
                    {
                        CRotPot(m, n, P, U, rmin, h, NumWFPoints);
                        RC = true;
                    }
                    fastWaveFunc(U, E[0][n][m][l], PM, h, NumWFPoints, WF, Ri = rmin, Ra = rmax, ISq, F, i,
                                 c, a);
                    F1 = F2 = 1.0 / (ISq * ISq);
                    F1 *= W1[n][l][m];
                    F2 *= W2[n][l][m];
                    for (k=i; k <= c; k++) 
                    {
                        WFS[k] += F1 * WF[k] * WF[k];
                        WWFS[k] += F2 * WF[k] * WF[k];
                    }
                    for (F1 *= F * F, F2 *= F * F; k<a; k++) 
                    {
                        WFS[k] += F1 * WF[k] * WF[k];
                        WWFS[k] += F2 * WF[k] * WF[k];
                    }
                }
            }
        }
        delete[] WF;
        delete[] P;
        delete[] U;
        Destroy(W1, NIso, Mv + 1);
        Destroy(W2, NIso, Mv + 1);
    }
    //printf("Vor 3.\n");
    for (k=0; k < NK; k++) if (KS[k] > KM)
    {
        KM = KS[k];
        *rv = k;
    }
    if (KS != 0) delete[] KS;
    int e, g;
    if (Ndev != N)
    {
        if (Ndev > 0)
        {
            delete[] Row;
            delete[] dev;
        }
        Row = new int[Ndev = N];
        dev = new double[N];
    }
    for (n=l=e=g=0; n<N; l++)
    {
        if (ES[l]) 
        {
            for (m=0; m < NEL[e]; m++)
            {
                Row[n] = EL[e][m].Row;
                dev[n++] = EL[e][m].dev;
                //DevR[n++] = EL[e][m].DevR;
            }
            e++;
        }
        else
        {
            for (m=0; m < nGL[g]; m++)
            {
                Row[n] = GL[g][m].Row;
                dev[n++] = GL[g][m].dev;
                //DevR[n++] = GL[g][m].DevR;
            }
            g++;
        }
    }
    BLLock->unlock();
    if (ShowBadList)
    {
        QString **Data;
        showBadList(Data, m, n, NumWFPoints);
        BLLock->unlock();
    }
    if (hCi_cP > 0.0)
    {
        bool *LRCFree = getLRCFree();
        double *LRC, Ra;
        int NLRC, *pLRC;
        getLRCoeff(Ra, NLRC, pLRC, LRC);
        for (n=0; n < NLRC; n++) if (LRCFree[n] && pLRC[n] >= 12)
        {
            CE = pow(LRC[n-1] / LRC[n-2], 3) * LRC[n-3];
            D = (CE - LRC[n]) / (hCi_cP * CE);
            FQS += D*D;
        }
    }
    for (n=0; n < NumWantedValues; n++)
    {
        D = (WantedValues[n].Value - WantedValues[n].CurValue) / WantedValues[n].Precision;
        FQS += D*D;
    }
    if (LevenbergMarquardtRunning) emit badListChanged(FQS, lambda, minLSFQS);
    else emit badListChanged(FQS, -1.0, 0.0);
    return FQS;
}

void PotWorker::getDiagFuncs(double*& rWFS, double*& rWWFS, int NumWFPoints)
{
    int n;
    if (WFS == 0)
    {
        rWFS = rWWFS = 0;
        return;
    }
    rWFS = new double[NumWFPoints];
    rWWFS = new double[NumWFPoints];
    BLLock->lock();
    for (n=0; n < NumWFPoints; n++)
    {
        rWFS[n] = WFS[n];
        rWWFS[n] = WWFS[n];
    }
    BLLock->unlock();
}

void PotWorker::showBadList(QString**& Data, int& NR, int &NC, int NumWFPoints)
{
    BLLock->lock();
    int NG, NE, s, n, m, i;
    if (dev == 0)
    {
        if(!Fit->isRunning()) 
        {
            int N = numSplinePoints + 2, **QN;
            int rv = N, F, I, J, NFC = int(2.0 * el_S) + 1;
            double ****E, *****K, **EQS, *EO, *DE, *Sig, ****MaxErr, *Err;
            estMinErr(MaxErr, NFC, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
            CreateEK(E, K, mJ, mv, N);
            CreateFitEQS(NE, EQS, EO, DE, Sig, Err, QN, N, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);
            calcE(E, mJ, mv, Mv, NumWFPoints, MaxErr);
            calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, false, 0, 0, K, &rv);
            for (F=0; F < NFC; F++)
            {
                for (I=0; I < NIso; I++)
                {
                    for (J=0; J <= mJ[I]; J++) if (MaxErr[F][I][J] != 0) delete[] MaxErr[F][I][J];
                    if (MaxErr[F][I] != 0) delete[] MaxErr[F][I];
                }
                delete[] MaxErr[F];
            }
            delete[] MaxErr;
            DestroyFitData(E, K, mJ, mv, NE, EQS, EO, DE, Sig, Err, QN);
        }
        else
        {
            Data = 0;
            NR = NC = 0;
            return;
        }
    }
    for (s = NG = 0; s < nGS; s++) for (n=0; n < nGL[s]; n++) 
            if (fabs(GL[s][n].DevR) > 3.0 || (GL[s][n].DevR == 0.0 && GL[s][n].dev != 0.0)) NG++;
    for (s = NE = 0; s < nEs; s++) for (n=0; n < NEL[s]; n++)
            if (fabs(EL[s][n].DevR) > 3.0 || EL[s][n].DevR == 0.0) NE++;
    TermEnergy **EP = 0, *EB = EL[0];
    TableLine **GP = 0, *GB = GL[0];
    //printf("Vor NE\n");
    if (NE > 0)
    {
        TermEnergy **EP1 = new TermEnergy*[NE];
        EP = new TermEnergy*[NE];
        for (s=m=0; s < nEs; s++) for (n=0; n < NEL[s]; n++) 
                if (fabs(EL[s][n].DevR) > 3.0 || EL[s][n].DevR == 0.0) EP1[m++] = EL[s] + n;
        //printf("NE=%d, NEL[0]=%d\n", NE, NEL[0]);
        while (EB != 0) for (n=0, m=2, EB = 0; m <= NE; n++, m+=2) 
        {
            //printf("n=%d, m=%d\n", n, m);
            if (m < NE ? ((EP1[m-1]->DevR == 0.0 && EP1[m]->DevR != 0.0) || (EP1[m-1]->DevR != 0.0 
                         && EP1[m]->DevR != 0.0 && fabs(EP1[m-1]->DevR) > fabs(EP1[m]->DevR))) : true)
                i=m-1;
            else i=m;
            //printf("i=%d\n", i);
            if ((EP1[i]->DevR == 0.0 && EP1[n]->DevR != 0.0) 
                         || (EP1[i]->DevR != 0.0 && EP1[n]->DevR != 0.0 
                         && fabs(EP1[i]->DevR) > fabs(EP1[n]->DevR)))
            {
                EB = EP1[i];
                EP1[i] = EP1[n];
                EP1[n] = EB;
            }
        }
        //printf("Nach S1\n");
        for (i=0; i < NE; i++)
        {
            EP[i] = EP1[0];
            for (n=0, m=2; EP1[n] != 0; m=2*((n=m)+1))
            {
                if (m < NE)
                {
                    if (EP1[m-1] != 0 ? (EP1[m] != 0 ? (EP1[m-1]->DevR == 0.0 && EP1[m]->DevR != 0.0) 
                        || (EP1[m-1]->DevR != 0.0 && EP1[m]->DevR != 0.0 
                        && fabs(EP1[m-1]->DevR) > fabs(EP1[m]->DevR)) : true) : false) m--;
                }
                else if (m==NE) m--;
                else break;
                EP1[n] = EP1[m];
            }
            if (m >= NE) EP1[n] = 0;
        }
        //printf("Nach S2\n");
        delete[] EP1;
        if (NG == 0)
        {
            Data = CreateQString(NE, 7);
            for (n=0; n < NE; n++)
            {
                Data[n][0] = QString::number(EP[n]->Iso);
                Data[n][1] = QString::number(EP[n]->v);
                Data[n][2] = QString::number(EP[n]->J);
                Data[n][3] = QString::number(EP[n]->E, 'f', EP[n]->nDig);
                Data[n][4] = QString::number(EP[n]->err, 'f', EP[n]->nDig);
                if (EP[n]->DevR != 0.0)
                {
                    Data[n][5] = QString::number(EP[n]->dev, 'f', EP[n]->nDig);
                    Data[n][6] = QString::number(EP[n]->DevR, 'f', 2);
                }
                else Data[n][5] = "not found";
            }
            BadListNC = NC = 7;
        }
    }
    //printf("Vor NG\n");
    if (NG > 0)
    {
        GP = new TableLine*[NG];
        TableLine **GP1 = new TableLine*[NG];
        for (s=m=0; s < nGS; s++) for (n=0; n < nGL[s]; n++) 
                if ((GL[s][n].DevR == 0.0 && GL[s][n].dev != 0.0) || fabs(GL[s][n].DevR) > 3.0) 
                  GP1[m++] = GL[s] + n;
        while(GB != 0) for (n=0, m=2, GB = 0; m <= NG; n++, m+=2) 
        {
            //printf("s=%d\n", s);
            if (m < NG ? (GP1[m-1]->DevR == 0.0 && GP1[m]->DevR != 0.0) || (GP1[m-1]->DevR != 0.0 
                         && GP1[m]->DevR != 0.0 && fabs(GP1[m-1]->DevR) > fabs(GP1[m]->DevR)) : true)
                i=m-1;
            else i=m;
            if ((GP1[i]->DevR == 0.0 && GP1[n]->DevR != 0.0) || 
                          (GP1[i]->DevR != 0.0 && GP1[n]->DevR != 0.0
                          && fabs(GP1[i]->DevR) > fabs(GP1[n]->DevR)))
            {
                //printf("GP[%d]=%f, GP[%d]=%f, NG=%d\n", i, GP1[i]->DevR, n, GP1[n]->DevR, NG);
                GB = GP1[i];
                GP1[i] = GP1[n];
                GP1[n] = GB;
            }
        }
        for (i=0; i < NG; i++)
        {
            //printf("i=%d, NG=%d\n", i, NG);
            GP[i] = GP1[0];
            for (n=0, m=2; GP1[n] != 0; m=2*((n=m)+1))
            {
                if (m < NG)
                {
                    if (GP1[m-1] != 0 ? (GP1[m] != 0 ? (GP1[m-1]->DevR == 0.0 && GP1[m]->DevR != 0.0) 
                        || (GP1[m-1]->DevR != 0.0 && GP1[m]->DevR != 0.0 
                                   && fabs(GP1[m-1]->DevR) > fabs(GP1[m]->DevR)) : true) : false) m--;
                }
                else if (m==NG) m--;
                else break;
                GP1[n] = GP1[m];
            }
            if (m >= NG) GP1[n] = 0;
        }
        delete[] GP1;
        //printf("Nach Sortieren\n");
        Data = CreateQString(NG + NE, 9);
        for (s=m=n=0; n < NG + NE; n++) 
        {
            //printf("m=%d, NE=%d, s=%d, NG=%d, n=%d\n", m, NE, s, NG, n);
            if (s < NG ? (m < NE ? ((fabs(GP[s]->DevR) > fabs(EP[m]->DevR) && EP[m]->DevR != 0.0)
                         || GP[s]->DevR == 0.0) : true) : false)
            {
                Data[n][0] = QString::number(GP[s]->Iso);
                Data[n][1] = QString::number(GP[s]->vss);
                Data[n][2] = QString::number(GP[s]->Jss);
                Data[n][3] = QString::number(GP[s]->vs);
                Data[n][4] = QString::number(GP[s]->Js);
                Data[n][5] = QString::number(GP[s]->WN, 'f', GP[s]->nDig);
                Data[n][6] = QString::number(GP[s]->err, 'f', GP[s]->nDig);
                if (GP[s]->DevR != 0.0)
                {
                    Data[n][7] = QString::number(GP[s]->dev, 'f', GP[s]->nDig);
                    Data[n][8] = QString::number(GP[s]->DevR, 'f', 2);
                }
                else Data[n][7] = "not found";
                s++;
            }
            else
            {
                Data[n][0] = QString::number(EP[m]->Iso);
                Data[n][1] = QString::number(EP[m]->v);
                Data[n][2] = QString::number(EP[m]->J);
                Data[n][5] = QString::number(EP[m]->E, 'f', EP[m]->nDig);
                Data[n][6] = QString::number(EP[m]->err, 'f', EP[m]->nDig);
                if (EP[m]->DevR != 0.0)
                {
                    Data[n][7] = QString::number(EP[m]->dev, 'f', EP[m]->nDig);
                    Data[n][8] = QString::number(EP[m]->DevR, 'f', 2);
                }
                else Data[n][7] = "not found";
                m++;
            }
        }
        BadListNC = NC = 9;
    }
    if (BadListNR > 0) Destroy(BadListData, BadListNR);
    BadListData = Data;
    BadListNR = NR = NE + NG;
    //printf("Nach NG\n");
    if (EP != 0) delete[] EP;
    if (GP != 0) delete[] GP;
}

void PotWorker::unlockBadList()
{
    BLLock->unlock();
}

void PotWorker::getBadList(QString**& Data, int& NR, int& NC)
{
    BLLock->lock();
    Data = BadListData;
    NR = BadListNR;
    NC = BadListNC;
}

void PotWorker::getDeviations(int& N, int*& rRow, double*& rdev, double*& rRWErr)
{
    BLLock->lock();
    N = Ndev;
    rRow = Row;
    rdev = dev;
    rRWErr = RWErr;
}

void PotWorker::badListChanged(QString* Items, int, TableLine &L)
{
    if (Fit->isRunning()) return;
    int n, m;
    if (nGS == 0 || Items[3].isEmpty()) 
    {
        for (n=0; n < nEs; n++) for (m=0; m < NEL[n]; m++) 
            if (EL[n][m].Iso == L.Iso && EL[n][m].v == L.vss && EL[n][m].J == L.Jss && EL[n][m].E == L.WN)
                EL[n][m].err = L.err;
    }
    else for (n=0; n < nGS; n++) for (m=0; m < nGL[n]; m++)
        if (GL[n][m].Iso == L.Iso && GL[n][m].vss == L.vss && GL[n][m].Jss == L.Jss && GL[n][m].vs == L.vs
                && GL[n][m].Js == L.Js && GL[n][m].WN == L.WN) GL[n][m].err = L.err;
}

void PotWorker::destroyFitData()
{
    //printf("destroyFitData\n");
    int I, n;
    for (I=0; I < NIso; I++) if (mJ[I] > -1) delete[] mv[I];
    //printf("Nach 0.\n");
    delete[] mv;
    delete[] mJ;
    for (n=0; n < nEs; n++) delete[] EL[n];
    //printf("Nach 1.\n");
    delete[] EL;
    delete[] NEL;
    for (n=0; n < nGS; n++) if (nGP[n] > 0)
    {
        delete[] PS[n];
        delete[] GL[n];
    }
    //printf("Nach zweiter Schleife\n");
    delete[] PS;
    delete[] GL;
    delete[] nGL;
    delete[] nGP;
    nEs = nGS = 0;
    //printf("Ende destroyFitData\n");
}

void PotWorker::DestroyFitData(double ****E, double *****K, int *mJ, int **mv, int NL, double **EQS, double *EO,
                               double *DE, double *Sig, double *Err, int *QN[5])
{
    //printf("Potential::DestroyFitData\n");
    int I, J, v, F, NFC = int(2.0 * el_S) + 1;
    for (F=0; F < NFC; F++)
    {
        for (I=0; I < NIso; I++) if (mJ[I] > -1)
        {
            //printf("I=%d, NIso=%d, mJ[I]=%d\n", I, NIso, mJ[I]);
            for (J=0; J <= mJ[I]; J++) if (mv[I][J] > -1)
            {
                //printf("J=%d, mJ[I]=%d, mv[I][J]=%d\n", J, mJ[I], mv[I][J]);
                delete[] E[F][I][J];
                for (v=0; v <= mv[I][J]; v++) delete[] K[F][I][J][v];
                delete[] K[F][I][J];
            }
            delete[] E[F][I];
            delete[] K[F][I];
        }
        delete[] E[F];
        delete[] K[F];
    }
    //printf("Nach erster Schleife\n");
    delete[] E;
    delete[] K;
    Destroy(EQS, NL + NumWantedValues);
    delete[] EO;
    delete[] DE;
    delete[] Sig;
    delete[] Err;
    Destroy(QN, NL);
    //printf("Ende\n");
}

void PotWorker::estMinErr(double****& MinErr, int NFC, int* mJ, int** mv, int nEs, int* NEL, TermEnergy** EL, 
                          int nGS, int* nGL, TableLine** GL)
{
    int F, I, J, v, s, l;
    double ErrF = 0.1, MaxErr = 0.001, EB;
    for (F = 0, MinErr = new double***[NFC]; F < NFC; F++) for (I=0, MinErr[F] = new double**[NIso]; I < NIso; I++)
    {
        if (mJ[I] >= 0) for (J=0, MinErr[F][I] = new double*[mJ[I] + 1]; J <= mJ[I]; J++) 
            for (v=0, MinErr[F][I][J] = new double[mv[I][J] + 1]; v <= mv[I][J]; v++) MinErr[F][I][J][v] = MaxErr;
        else MinErr[F][I] = 0;
    }
    for (s=0; s < nGS; s++) for (l=0; l < nGL[s]; l++)
    {
        F = (GL[s][l].FC < NFC && GL[s][l].FC >= 0 ? GL[s][l].FC : 0);
        if ((EB = ErrF * GL[s][l].err) < MinErr[F][GL[s][l].Iso][GL[s][l].Jss][GL[s][l].vss])
            MinErr[F][GL[s][l].Iso][GL[s][l].Jss][GL[s][l].vss] = EB;
    }
    for (s=0; s < nEs; s++) for (l=0; l < NEL[s]; l++)
    {
        F = (EL[s][l].FC < NFC && EL[s][l].FC >= 0 ? EL[s][l].FC : 0);
        if ((EB = ErrF * EL[s][l].err) < MinErr[F][EL[s][l].Iso][EL[s][l].J][EL[s][l].v])
            MinErr[F][EL[s][l].Iso][EL[s][l].J][EL[s][l].v] = EB;
    }
}

double PotWorker::getFQS(int NumWFPoints, bool calc, double *StdDev, bool Thread, double SFQSRad, double *sigma, double *FQS_PAL)
{
    if (mJ == 0) return -1.0;
    int n, N, m, I, J, NFC = int(el_S * 2.0) + 1, F, MJ = 0;
    double ****E = new double***[NFC], ****MaxErr;
    if (SFQSRad > 0.0) SFQSU = new bool***[NFC];
    if (!Fit->isRunning() || Thread)
    {
        LevenbergMarquardtRunning = false;
        for (F = 0; F < NFC; F++)
        {
            E[F] = new double**[NIso];
            if (SFQSRad > 0.0) SFQSU[F] = new bool**[NIso];
            for (I=0; I < NIso; I++) 
            {
                if (mJ[I] > -1)
                {
                    E[F][I] = new double*[mJ[I] + 1];
                    if (SFQSRad > 0.0) SFQSU[F][I] = new bool*[mJ[I] + 1];
                    for (J=0; J <= mJ[I]; J++)
                    {
                        E[F][I][J] = new double[mv[I][J] + 1];
                        if (SFQSRad > 0.0) SFQSU[F][I][J] = new bool[mv[I][J] + 1];
                    }
                }
                else E[F][I] = 0;
            }
        }
        //printf("Vor calcE\n");
        estMinErr(MaxErr, NFC, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
        calcE(E, mJ, mv, Mv, NumWFPoints, MaxErr, (SFQSRad > 0.0 ? SFQSU : 0), SFQSRad);
        //printf("Vor calcFQS\n");
        if (calc) for (n=0; n < NIso; n++) if (mJ[n] > MJ) MJ = mJ[n];
        FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, calc, MJ, Mv, 0, 0, 0, 1.0, 0.0, SFQSU);
        //printf("Vor Delete\n");
        if (StdDev != 0)
        {
            for (n=N=0, *StdDev = 0.0; n < nEs; n++) for (m=0; m < NEL[n]; m++, N++) *StdDev += EL[n][m].dev * EL[n][m].dev;
            for (n=0; n < nGS; n++) for (m=0; m < nGL[n]; m++, N++) *StdDev += GL[n][m].dev * GL[n][m].dev;
            *StdDev = sqrt(*StdDev / (N - (NumFreePar != 0 ? NumFreePar : 1)));
        }
        if (FQS_PAL != 0)
        {
            for (n=0, *FQS_PAL = 0.0; n < nEs; n++) for (m=0; m < NEL[n]; m++) if (EL[n][m].err < 1e-5) 
                *FQS_PAL += EL[n][m].DevR * EL[n][m].DevR;
            for (n=0; n < nGS; n++) for (m=0; m < nGL[n]; m++) if (GL[n][m].err < 1e-5) *FQS_PAL += GL[n][m].DevR * GL[n][m].DevR;
        }
        for (F=0; F < NFC; F++)
        {
            for (I=0; I < NIso; I++)
            {
                //printf("I=%d, NIso=%d\n", I, NIso);
                for (J=0; J<= mJ[I]; J++) 
                {
                    if (E[F][I][J] != 0) delete[] E[F][I][J];
                    if (MaxErr[F][I][J] != 0) delete[] MaxErr[F][I][J];
                }
                if (E[F][I] != 0) delete[] E[F][I];
                if (MaxErr[F][I] != 0) delete[] MaxErr[F][I];
                if (mv[I] != 0 && SFQSRad <= 0.0) 
                {
                    delete[] mv[I];
                    mv[I] = 0;
                }
            }
            delete[] E[F];
            delete[] MaxErr[F];
        }
        delete[] E;
        delete[] MaxErr;
        if (SFQSRad <= 0.0)
        {
            delete[] mJ;
            delete[] mv;
        }
        for (n=0; n < nEs; n++) delete[] EL[n];
        delete[] EL;
        delete[] NEL;
        for (n=0; n < nGS; n++) if (nGL[n] > 0)
        {
            delete[] PS[n];
            delete[] GL[n];
        }
        delete[] PS;
        delete[] GL;
        delete[] nGL;
        delete[] nGP;
        nEs = nGS = 0;
    }
    else if (StdDev != 0)
    {
        BLLock->lock();
        for (n=0, *StdDev = 0.0; n < Ndev; n++) *StdDev += dev[n] * dev[n];
        *StdDev = sqrt(*StdDev / (Ndev - NumFreePar));
        BLLock->unlock();
    }
    if (sigma != 0) *sigma = sqrt(FQS / (Ndev - NumFreePar));
    return FQS;
}

void PotWorker::getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients)
{
    R = (numSplinePoints > 0 ? points[numSplinePoints - 1].x : 0.0);
    numCoefficients = NLRC;
    Exponents = PLRC;
    Coefficients = LRC;
}

void PotWorker::setLRCoeff(double R, int numCoefficients, int* Exponents, double* Coefficients, bool* LRCFree)
{
    if (R > 0.0 && numSplinePoints > 0 && R < points[numSplinePoints - 1].x) Ro = R;
    NLRC = numCoefficients;
    if (PLRC != 0 && PLRC != Exponents) delete[] PLRC;
    PLRC = Exponents;
    if (LRC != 0 && LRC != Coefficients) delete[] LRC;
    LRC = Coefficients;
    if (LRCfree != 0 && LRCfree != LRCFree) delete[] LRCfree;
    LRCfree = LRCFree;
}

void PotWorker::createSMap()
{
    int n, m;
    for (NumFreePar = n = 0; n < numSplinePoints; n++) if (points[n].variable) NumFreePar++;
    for (NumFreeLRC = n = 0; n < NLRC; n++) if (LRCfree[n]) NumFreeLRC++;
    //if (State->getStateNum() > 0 && NumFreeLRC > 0) NumFreeLRC = 1;
    for (NumFreeAdCorr = n = 0; n < NAdCorr; n++) if (adCorrFree[n]) NumFreeAdCorr++;
    NumFreePar += (NumFreeLRC + NumFreeAdCorr);
    if (shortRVariable) NumFreePar++;
    SMap = new int[NumFreePar];
    if (shortRVariable) 
    {    
        SMap[0] = 0;
        n=1;
    }
    else n=0;
    for (m=0; m < numSplinePoints; m++) if (points[m].variable) SMap[n++] = m + 1;
    if (NumFreeLRC > 0) 
    {
        /*if (State->getStateNum() > 0) SMap[n++] = numSplinePoints + 1;
        else*/ for (m=0; m < NLRC; m++) if (LRCfree[m]) SMap[n++] = numSplinePoints + 1 + m;
    }
    if (NumFreeAdCorr > 0) for (m=0; m < NAdCorr; m++) if (adCorrFree[m]) 
        SMap[n++] = numSplinePoints + NLRC + m + 1;
}

void PotWorker::saveCoefficients(double*& bC)
{
    int n, m;
    if (bC == 0) bC = new double[2 * numSplinePoints + NLRC + NAdCorr + 2];
    bC[0] = iA;
    bC[1] = iO;
    for (n=0, m=2; n < numSplinePoints; n++)
    {
        bC[m++] = points[n].y;
        bC[m++] = points[n].yss;
    }
    for (n=0; n < NLRC; n++) bC[m++] = LRC[n];
    for (n=0; n < NAdCorr; n++) bC[m++] = adCorr[n];
}

void PotWorker::restorePotential(double* bC)
{
    int n, m;
    Lock->lock();
    iA = bC[0];
    iO = bC[1];
    for (n=0, m=2; n < numSplinePoints; n++)
    {
        points[n].y = bC[m++];
        points[n].yss = bC[m++];
    }
    for (n=0; n < NLRC; n++) LRC[n] = bC[m++];
    for (n=0; n < NAdCorr; n++) adCorr[n] = bC[m++];
    Lock->unlock();
}

void PotWorker::setCurValuesForWantedValues()
{
    int n;
    for (n=0; n < NumWantedValues; n++)
    {
        if (WantedValues[n].CoeffNum == 0) WantedValues[n].CurValue = -iExp * iA * pow(points[0].x, -iExp - 1);
        else if (WantedValues[n].CoeffNum <= numSplinePoints) WantedValues[n].CurValue = points[WantedValues[n].CoeffNum - 1].y;
        else if (WantedValues[n].CoeffNum <= numSplinePoints + NLRC) WantedValues[n].CurValue = LRC[WantedValues[n].CoeffNum - numSplinePoints - 1];
        else WantedValues[n].CurValue = adCorr[WantedValues[n].CoeffNum - numSplinePoints - NLRC - 1];
    }
}

void PotWorker::FitSplinePot(int NumWFPoints, int MaxIt, double threshFakt)
{
    int n, r, NFC, F, I, J;
    double FQSv = -1.0, FQSn, minR, maxR, ****MaxErr;
    bool RPoints[numSplinePoints];
    LevenbergMarquardtRunning = false;
    emit fitRunning();
    if (SMap == 0) createSMap();
    for (n=0; n < NumFreePar; n++) printf("SMap[%d]=%d\n", n, SMap[n]);
    estMinErr(MaxErr, NFC = int(2.0 * el_S) + 1, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
    //for (n=0; n < nEs; n++) NLines += NEL[n];
    //for (n=0; n < nGS; n++) NLines += nGL[n];

    while ((r = splineFit(RPoints, FQSv, FQSn, mJ, mv, Mv, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, 
            minR, maxR, NumWFPoints, MaxErr, MaxIt, threshFakt)) != -1)
    {
        break;
        
        /*if (FQSn > 0.999 * lFQS) break;
        lFQS = FQSn;
        printf("Optimize points\n");
        r=-2;
        for (n=N=0; n < numSplinePoints; n++) if (RPoints[n]) N++;
        if (r>=0) N++;
        PB = new SplinePoint[N];
        for (n=m=0; n < numSplinePoints; n++)
        {
            if (n==r)
            {
                switch (r)
                {
                    case 0:
                        PB[0].x = ((X = 2.0 * points[0].x - points[1].x) > minR ? X : minR);
                        PB[m++].y = getPoint(PB[0].x);
                        break;
                    case 1:
                        PB[0].x = ((X = 1.5 * points[0].x - 0.5 * points[1].x) > minR ? X : minR);
                        PB[0].y = getPoint(PB[0].x);
                        PB[1].x = PB[0].x + 0.75 * (points[0].x - X);
                        PB[m++].y = getPoint(PB[1].x);
                        break;
                    default:
                        PB[m-1].x = 0.25 * points[n-2].x + 0.75 * points[n-1].x;
                        PB[m-1].y = getPoint(PB[m-1].x);
                        PB[m].x = 0.75 * points[n-1].x + 0.25 * points[n].x;
                        PB[m].y = getPoint(PB[m].x);
                        m++;
                        break;
                }
            }
            if (RPoints[n]) PB[m++] = points[n];
        }
        //printf("m=%d, N=%d\n", m, N);
        if (r==n)
        {
            PB[N-1].x = ((X = 1.5 * points[numSplinePoints - 1].x 
                            - 0.5 * points[numSplinePoints - 2].x) < maxR ? X : maxR);
            PB[N-1].y = getPoint(PB[N-1].x);
            PB[N-2].x = PB[N-1].x + 0.75 * (points[numSplinePoints - 1].x - X);
            PB[N-2].y = getPoint(PB[N-2].x);
        }
        else if (r==n+1)
        {
            PB[N-1].x = ((X = 2.0 * points[numSplinePoints - 1].x - points[numSplinePoints - 2].x)
                          < maxR ? X : maxR);
            PB[N-1].y = getPoint(PB[N-1].x);
        }
        Ch = true;
        //printf("Vor zweiter\n");
        while (Ch) for (Ch = false, n=1; n<N-1; n++) 
        {
            if (PB[n].x - PB[n-1].x < 0.5 * (PB[n+1].x - PB[n].x))
            {
                PB[n].x = 0.66 * PB[n-1].x + 0.34 * PB[n+1].x;
                PB[n].y = getPoint(PB[n].x);
                Ch = true;
            }
            else if (PB[n].x - PB[n-1].x > 2.0 * (PB[n+1].x - PB[n].x))
            {
                PB[n].x = 0.34 * PB[n-1].x + 0.66 * PB[n+1].x;
                PB[n].y = getPoint(PB[n].x);
                Ch = true;
            }
        }
        if (N==0)
        {
            printf("Error: removed all points!\n");
            delete[] PB;
            break;
        }
        //printf("Vor delete, numSplinePoints=%d\n", numSplinePoints);
        Destroy(L, numSplinePoints);
        //printf("Vor delete S, K=%d\n", numSplinePoints+2);
        Destroy(S, numSplinePoints + 2);
        L = 0;
        S = 0;
        iA = 0.0;
        delete[] points;
        points = PB;
        numSplinePoints = N;
        delete[] LRC;
        delete[] PLRC;
        NLRC = 0;
        //printf("Vor calcyss\n");
        calcyss();
        emit propertiesChanged();
        printf("Ende optimization\n");
        if (!DoubleWell)
        {
            getMinimum(MR, MV);
            for (C = false, n=0; (n < numSplinePoints ? points[n].x < MR : false); n++) ;
            for (LV = MV; n < numSplinePoints; n++)
            {
                if ((AV = points[n].y) < LV)
                {
                    for (m=n-1; points[m].y > AV; m--) points[m+1].y = points[m].y;
                    points[m+1].y = AV;
                    C = true;
                }
                else LV = AV;
            }
            for (n--; (n >= 0 ? points[n].x > MR : false); n--) ;
            for (LV = MV; n >= 0; n--)
            {
                if ((AV = points[n].y) < LV)
                {
                    for (m=n+1; points[m].y > AV; m++) points[m-1].y = points[m].y;
                    points[m-1].y = AV;
                    C = true;
                }
                else LV = AV;
            }
            if (C) 
            {
                calcyss();
                emit propertiesChanged();
            }
        }*/
    }
    for (F=0; F < NFC; F++)
    {
        for (I=0; I < NIso; I++)
        {
            for (J=0; J <= mJ[I]; J++) if (MaxErr[F][I][J] != 0) delete[] MaxErr[F][I][J];
            if (MaxErr[F][I] != 0) delete[] MaxErr[F][I];
        }
        delete[] MaxErr[F];
    }
    delete[] MaxErr;
}

void PotWorker::updatePotential(double* C)
{
    Lock->lock();
    int i;
    double Rb, ys, Ra = points[numSplinePoints - 1].x;
    for (i=0; i < NumFreePar; i++) if (SMap[i] > 0 && SMap[i] <= numSplinePoints) 
        points[SMap[i] - 1].y += C[i];
    if (SMap[0] == 0) iA += (C[0] * pow(points[0].x, iExp));
    /*}
    else
    {
        double Par[numSplinePoints + NLRC + 1];
        for (n=0, Par[0] = iA; n < numSplinePoints; n++) Par[n+1] = points[n].y;
        for (n=0; n < NLRC; n++) Par[numSplinePoints + n + 1] = LRC[n];
        FQSn = fitSplinePot(this, E, K, mJ, mv, Mv, N, NE, QN, Sig, EQS, EO, DE, minR, maxR, 
                            numSplinePoints + NLRC + 1, Par);
        for (iA = Par[n=0]; n < numSplinePoints; n++) points[n].y = Par[n+1];
        for (n=0; n < NLRC; n++) LRC[n] = Par[numSplinePoints + n + 1];
    }*/
    /*if (GS) 
    {
        for (n=0, i=1; i < numSplinePoints; i++) if (points[i].y < points[n].y) n=i;
        if (points[numSplinePoints - 1].x > maxR) Uinf -= points[n].y;
        for (i=0; i < numSplinePoints; i++) points[i].y -= points[n].y;
    }*/
    iO = points[0].y - iA * pow(points[0].x, -iExp);
    if (NLRC == 0)
    {
        LRC = new double[NLRC = 2];
        PLRC = new int[NLRC];
        PLRC[0] = 6;
        PLRC[1] = 8;
        LRC[0] = 0.0;
        LRC[1] = 0.0;
        Rb = pow(Ra, PLRC[0]);
        ys = PLRC[0] * LRC[0] * pow(Ra, -1 - PLRC[0]) 
           + PLRC[1] * LRC[1] * pow(Ra, -1 - PLRC[1]) 
           + (SMap[NumFreePar - 1] > numSplinePoints ? C[NumFreePar - 1] : 0.0);
        LRC[0] = Rb * (0.5 * double(2 + PLRC[0]) * (Uinf - points[numSplinePoints - 1].y) 
                       - 0.5 * (Ra) * ys);
        Rb *= Ra * Ra;
        LRC[1] = Rb * (0.5 * double(PLRC[0]) * (points[numSplinePoints - 1].y - Uinf) 
                + 0.5 * Ra * ys);
    }
    else 
    {
        for (i=0; i < NumFreePar; i++)
        {
            if (SMap[i] > numSplinePoints + NLRC) adCorr[SMap[i] - numSplinePoints - NLRC - 1] += C[i];
            else if (SMap[i] > numSplinePoints) 
                LRC[SMap[i] - numSplinePoints - 1] += (C[i] * pow(Ra, PLRC[SMap[i] - numSplinePoints - 1]));
        }
    }
    calcyss();
    Lock->unlock();
}

void PotWorker::calcyss()
{
    double ED, Rb, ys, Ra = points[numSplinePoints - 1].x;
    int i, n;
    if (L == 0) calcLMatrix();
    if (Ro == 0.0)
    {
        for (i=0, ED = points[numSplinePoints - 1].y, ys = 0.0; i < NLRC; i++)
        {
            ED += (Rb = LRC[i] * pow(Ra, -PLRC[i]));
            ys += double(PLRC[i]) / Ra * Rb;
        }
        ED -= Uinf;
        //if (GS)
        //{
        for (n=0, iO -= ED; n < numSplinePoints; n++) points[n].y -= ED;
        //Uinf = 0.0;
        //}
        for (n=0; n < numSplinePoints; n++) 
        {
            for (i=0, points[n].yss = iA * L[n][0]; i < numSplinePoints; i++) 
                points[n].yss += points[i].y * L[n][i+1];
            points[n].yss += ys * L[n][numSplinePoints + 1];
            //printf("YSS[%d]=%f\n", n, F);
        }
    }
    else for (n=0; n < numSplinePoints; n++) 
        for (i=0, points[n].yss = 0.0; i < numSplinePoints; i++)
            points[n].yss += points[i].y * L[n][i];
}

int PotWorker::splineFit(bool *RPoints, double &FQSv, double &FQSn, int *&mJ, int **&mv, int &Mv, 
                         int &nEs, int *&NEL, TermEnergy **&EL, int &nGS, int *&nGP, int *&nGL, 
                         int **&PS, TableLine **&GL, double &minR, double &maxR, int NumWFPoints, double ****MaxErr, int MaxIt, double threshFakt)
{
    //printf("Potential::splineFit\n");
    int i, it, n, **QN;
    int N = NumFreePar, NE, BadC;
    int rv = N, MJ;
    for (n = MJ = 0; n < NIso; n++) if (mJ[n] > MJ) MJ = mJ[n];
    double ****E, *****K, **EQS, *EO, *DE, *Sig, C[N], aFQS = 0.0, *Err;
    FQS = FQSv = -1.0;
    //bool GS = (Uinf == 0);
    //for (n=0; n < nEs; n++) if (NEL[n] > 0) GS = false;
    CreateEK(E, K, mJ, mv, N);
    CreateFitEQS(NE, EQS, EO, DE, Sig, Err, QN, N, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);
    if (Type == SplinePotential && L==0) calcLMatrix();
    else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
    if (S==0) calcS(NumWFPoints, rmin, rmax);
    for (it = BadC = 0; it < MaxIt && BadC < (FQS < FQSv ? 3 : 5); it++)
    {
        aFQS = FQS;
        Lock->lock();
        calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr,
               RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, 
               NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2, mJ, mv, Mv, N, NE, QN, Sig, EQS, EO, DE, minR, maxR, NumWFPoints, MaxErr);
        Lock->unlock();
        //printf("Nach calcEK\n");
        FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints);
        if (FQSv == -1.0) FQSv = FQS;
        if (FQS < threshFakt * aFQS) BadC = 0;
        else BadC++;
        //if (NE < 20000)
        //{
        SvdFit(EQS, C, DE, Sig, N, NE, 1e-2);//-2
            //printf("Nach SvdFit\n");
            /*
            for (i=0; points[i].x < minR; i++) ;
            for (i--; i >= 0; i--)
            {
                if (C[i+1] * C[i+2] <= 0.0) C[i+1] = 0.0;
                else if (fabs(C[i+1] * points[i+1].y) > fabs(C[i+2] * points[i].y))
                C[i+1] = C[i+2] * points[i].y / points[i+1].y;
            }
            for (i = numSplinePoints - 1; points[i].x > maxR; i--) ;
            for (i++; i < numSplinePoints; i++)
            {
                if (C[i+1] * C[i] <= 0.0) C[i+1] = 0.0;
                else if (fabs(C[i+1] * (Uinf - points[i-1].y)) > fabs(C[i] * (Uinf - points[i].y)))
                    C[i+1] = C[i] * (Uinf - points[i].y) / (Uinf - points[i-1].y);
            }*/
        updatePotential(C);
        //emit propertiesChanged();
        /*if ((D1 = getPoint(D2 = 0.99 * points[numSplinePoints - 1].x 
                   + 0.01 * points[numSplinePoints - 2].x))
                  > (D3 = points[numSplinePoints - 1].y)) printf("RemP\n");*/
        //printf("x1=%f, y1=%f, x2=%f, y2=%f\n", D2, D1, points[numSplinePoints-1].x, D3);
    }
    //printf("Nach Fitschleife, minR=%f, maxR=%f\n", minR, maxR);
    //if (NE < 20000)
    //{
    calcE(E, mJ, mv, Mv, NumWFPoints, MaxErr);
    FQS = FQSn = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, true, MJ, Mv, K, &rv);
    //}

    if (FQSn < NE) rv = -1;
    for (i=0; (i < numSplinePoints ? points[i].x < minR : false); i++) RPoints[i] = false;
    while (i < numSplinePoints ? points[i].x <= maxR : false) RPoints[i++] = true;
    /*if (FQS < 0.9999 * aFQS) 
    {
        if (i == numSplinePoints) i--;
        rv = -2;
    }
    printf("i=%d, numSplinePoints=%d\n", i, numSplinePoints);
    while (i < numSplinePoints) RPoints[i++] = false;*/
    DestroyFitData(E, K, mJ, mv, NE, EQS, EO, DE, Sig, Err, QN);
    if (rv == -1) for (n=0; n < numSplinePoints; n++) if (!RPoints[n]) return -2;
    return rv;
}

void PotWorker::FitSplinePotToFunction(double* E, int NPoints, double MinR, double MaxR)
{
    int n, m, N;
    double R, h = (MaxR - MinR) / (NPoints - 1), *D = new double[NPoints], **EQS, aFQS = 1e99 , FQS = 1e98, *C, *Sig = new double[NPoints];
    if (L == 0) calcLMatrix();
    calcS(NPoints, MinR, MaxR);
    if (SMap == 0) createSMap();
    C = new double[NumFreePar];
    for (N=0; (N < NumFreePar ? SMap[N] <= numSplinePoints + NLRC : false); N++) ;
    EQS = Create(NPoints, N);
    for (n=0; n < NPoints; n++) for (m=0, Sig[n] = 1e4; m < N; m++) EQS[n][m] = S[SMap[m]][n];
    while (FQS < aFQS)
    {
        aFQS = FQS;
        for (n=0, FQS = 0.0, R = MinR; n < NPoints; n++, R += h)
        {
            D[n] = E[n] - Point(R);
            FQS += D[n] * D[n];
        }
        printf("FSPtF: FQS=%f\n", FQS * 1E8);
        SvdFit(EQS, C, D, Sig, N, NPoints, 1e-2);
        for (n=N; n < NumFreePar; n++) C[n] = 0.0;
        updatePotential(C);
    }
    delete[] D;
    Destroy(EQS, NPoints);
    Destroy(S, NumS);
    S=0;
    delete[] Sig;
    delete[] C;
}

void PotWorker::SetSplinePar(double* Par, double ****E, double *****K, int *mJ, int **mv, int Mv, int N, 
                             int NL, int **QN, double *Sig, double **EQS, double *EO, double *DE, 
                             double &minR, double &maxR, int NumWFPoints)
{
    int n, m;
    double ys, Rb;
    iA = Par[0];
    iO = Par[1] - pow(points[0].x, -6.0);
    for (n=0; n < numSplinePoints; n++) points[n].y = Par[n+1];
    for (n=0, Uinf = ys = 0.0; n < NLRC; n++)
    {
        Uinf += (Rb = (LRC[n] = Par[numSplinePoints + n + 1]) * pow(points[numSplinePoints - 1].x, -PLRC[n]));
        ys += double(PLRC[n]) * Rb / points[numSplinePoints - 1].x;
    }
    for (n=0; n < numSplinePoints; n++)
    {
        for (m=0, points[n].yss = iA * L[n][0]; m < numSplinePoints; m++)
            points[n].yss += points[m].y  * L[n][m+1];
        points[n].yss += ys * L[n][m+1];
    }
    if (Type == SplinePotential && L==0) calcLMatrix();
    else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
    if (S==0) calcS(NumWFPoints, rmin, rmax);
    Lock->lock();
    calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr, RIso1AdCorr,
           RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, NSpinRGamma, SpinRGamma,
           SpinRR1, SpinRR2, mJ, mv, Mv, N, NL, QN, Sig, EQS, EO, DE, minR, maxR, NumWFPoints);
    Lock->unlock();
}

void PotWorker::MonteCarloSimIteration(double* MCSTab, double UncFact)
{
    int e, n;
    srand(time(0));
    for (e=0; e < nEs; e++) for (n=0; n < NEL[e]; n++) EL[e][n].E -= EL[e][n].dev + MCSTab[rand() % 32768] * EL[e][n].err * UncFact;
    for (e=0; e < nGS; e++) for (n=0; n < nGL[e]; n++) GL[e][n].WN -= GL[e][n].dev + MCSTab[rand() % 32768] * GL[e][n].err * UncFact;
    /*if (points != 0) FitSplinePot(100);
    else if (Type == analyticalPotential)*/ 
    FitAnaPot(false, true, (Type == SplinePotential ? 0.0 : 1.0), false);
}

void PotWorker::getFitResult(double& rFQS, double& rWeightSum, int& rNWeightF, int*& rWeightFRow, double*& rWeightF, int &numFreePar)
{
    BLLock->lock();
    rFQS = FQS;
    rWeightSum = WeightSum;
    rWeightF = WeightF;
    rNWeightF = NWeightF;
    rWeightFRow = WeightFRow;
    numFreePar = NumFreePar;
}

int PotWorker::cAdRows()
{
    int n, AdRows;
    for (n = AdRows = 0; n < NLRC; n++) if (LRCfree[n] && PLRC[n] >= 12) AdRows++;
    return AdRows;
}

void PotWorker::fillWantedValueEQS(int NEQ, double* Sig, double** EQS, double* DE)
{
    int n, m, c;
    setCurValuesForWantedValues();
    for (m = NEQ, n=0; n < NumWantedValues; n++, m++)
    {
        DE[m] = WantedValues[n].Value - WantedValues[n].CurValue;
        Sig[m] = 1.0 / WantedValues[n].Precision;
        //if (Sig[m] * DE[m] > 100.0) Sig[m] = 100.0 / DE[m];
        for (c=0; c < NumFreePar; c++) EQS[m][c] = (SMap[c] == WantedValues[n].CoeffNum ? 1.0 : 0.0);
    }
}

double PotWorker::FitAnaPot(int NumWFPoints, bool robustWeighting, bool UseSVD, double threshhold, bool useLevenbergMarquardt, int MaxIt,
                            bool adjustTEWeightFact, double hCi_cP, bool FitRestarted)
{
    int n, m, e, NLevels = 0, NEQ, **QN, F, I, J, BadCount, i, it=1, done, AdRows;
    double ****MaxErr, ****E, *****K, **EQS, *EObs, *EDiff, *Sig, *Err;
    double aFQS =9e99, minR, maxR, FQSv = -1.0, **LMEQS, **tLMEQS, NFC;

    //writePoints("DebugPointsStart.Dat", 3.0, 50.0, 4701, 8, false, true, 0);
    
    StopFit = FitStopped = false;
    if (!FitRestarted) 
    {
        LevenbergMarquardtRunning = false;
        FQS = 8e99;
    }
    emit fitRunning();
    minLSFQS = -1.0;
    if (SMap != 0) delete[] SMap;
    createSMap();
    double C[NumFreePar];// aimp;
    double *bC = 0;
    cdConnectLR1C();
    cdConnectSR();
    if (getType() == SplinePotential) calcyss();
    saveCoefficients(bC);
    estMinErr(MaxErr, NFC = int(2.0 * el_S) + 1, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
    for (n=0; n < nEs; n++) NLevels += NEL[n];
    if (WeightF != 0) delete[] WeightF;
    WeightF = (adjustTEWeightFact ? new double[NLevels] : 0);
    if (nEs + nGS == 0) return 0.0; 
    if (hCi_cP > 0.0) AdRows = cAdRows();
    else AdRows = 0;
    CreateEK(E, K, mJ, mv, NumFreePar);
    CreateFitEQS(NEQ, EQS, EObs, EDiff, Sig, Err, QN, NumFreePar, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, AdRows);
    for (e=m=0; e < nEs; e++) for (n=0; n < NEL[e]; n++) Sig[m++] *= sqrt(EL[e][n].WeightFactor);
    if (Type == SplinePotential && L==0) calcLMatrix();
    else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
    if (S==0) calcS(NumWFPoints, rmin, rmax);
    Lock->lock();
    calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr, RIso1AdCorr,
           RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, NSpinRGamma, SpinRGamma,
           SpinRR1, SpinRR2, mJ, mv, Mv, NumFreePar, NEQ, QN, Sig, EQS, EObs, EDiff, minR, maxR, NumWFPoints, MaxErr);
    Lock->unlock();
    if (NumWantedValues > 0) fillWantedValueEQS(NEQ, Sig, EQS, EDiff);
    FQSv = aFQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, false, 0, 0, 0, 0, (robustWeighting ? Sig : 0), 1.0, hCi_cP);
    calcFitResult();
    emit firstFCFCalculated(Result, FQSv);
    
    //writePoints("DebugPointsStart2.Dat", 3.0, 50.0, 4701, 8, false, true, 0);
    
    if (UseSVD && !(LevenbergMarquardtRunning && FitRestarted)) for (BadCount =0; BadCount < 3; it++)  
    {
        SvdFit(EQS, C, EDiff, Sig, NumFreePar, NEQ + NumWantedValues, 1e-2);
        updatePotential(C);
        if (it == MaxIt || isStopFit())
        {
            if (it == 1) FQS = FQSv;
            useLevenbergMarquardt = false;
            break;
        }
        if (Type == SplinePotential && L==0) calcLMatrix();
        else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
        if (S==0) calcS(NumWFPoints, rmin, rmax);
        Lock->lock();
        calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr,
               RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, NSpinRGamma,
               SpinRGamma, SpinRR1, SpinRR2, mJ, mv, Mv, NumFreePar, NEQ, QN, Sig, EQS, EObs, EDiff, minR, maxR, NumWFPoints, MaxErr);
        Lock->unlock();
        if (NumWantedValues > 0) fillWantedValueEQS(NEQ, Sig, EQS, EDiff);
        FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, false, 0, 0, 0, 0, (robustWeighting ? Sig : 0), 1.0, hCi_cP);
        if (FQS < aFQS && FQS > 0.0 && !isnan(FQS))
        {
            if (FQS < 0.9 * aFQS || (!useLevenbergMarquardt && aFQS - FQS > threshhold)) BadCount = 0;
            else BadCount++;
            aFQS = FQS;
            saveCoefficients(bC);
        }
        else
        {
            BadCount++;
            //if (stopOnFirstStepFailure && FQS > FQSv) return 0.0;
        }
    }
    else FQS = aFQS;
    if (FQS > aFQS || FQS == 0.0 || isnan(FQS))
    {
        restorePotential(bC);
        if (Type == SplinePotential && L==0) calcLMatrix();
        else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
        if (S==0) calcS(NumWFPoints, rmin, rmax);
        Lock->lock();
        calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr,
               RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, NSpinRGamma,
               SpinRGamma, SpinRR1, SpinRR2, mJ, mv, Mv, NumFreePar, NEQ, QN, Sig, EQS, EObs, EDiff, minR, maxR, NumWFPoints, MaxErr);
        Lock->unlock();
        if (NumWantedValues > 0) 
        {
            for (n=0; n < NumFreePar; n++) C[n] = bC[SMap[n]];
            fillWantedValueEQS(NEQ, Sig, EQS, EDiff);
        }
    }
    if (useLevenbergMarquardt && !isStopFit())
    {
        printf("Now Levenberg-Marquardt:\n");  // Using knowledge from the book "Numerical Recipes, Third Edition"
        LMEQS = Create(NumFreePar, NumFreePar + 1);
        tLMEQS = Create(NumFreePar, NumFreePar + 1);
        calcLMEQS(NEQ + NumWantedValues, EQS, Sig, EDiff, LMEQS, hCi_cP);
        if (!(FitRestarted && LevenbergMarquardtRunning)) lambda = 1e-3;
        LevenbergMarquardtRunning = true;
        for (done = BadCount = 0; true; it++)
        {
            if (done == 4)
            {
                if (adjustTEWeightFact) for (m=e=0; e < nEs; e++) for (n=0; n < NEL[e]; n++) 
                    WeightF[m++] = (fabs(EL[e][n].DevR) < 4.0 ? (fabs(EL[e][n].DevR) < 1.0 ? 0.5 : 1.0) : 2.0);
                lambda = 0.0;
                
                //writePoints("DebugPointsFin.Dat", 3.0, 50.0, 4701, 8, false, true, 0);
            }
            if (done < 5)
            {
                for (n=0; n < NumFreePar; n++)
                {
                    for (i=0; i <= NumFreePar; i++) tLMEQS[n][i] = LMEQS[n][i];
                    tLMEQS[n][n] *= (1.0 + lambda);
                }
                SolvLinEqSwItImp(tLMEQS, NumFreePar, &minLSFQS);
                for (n=0; n < NumFreePar; n++) C[n] = tLMEQS[n][NumFreePar];
                updatePotential(C);
            }
            if (it == MaxIt || isStopFit())
            {
                if (it == 1) FQS = FQSv;
                break;
            }
            if (Type == SplinePotential && L==0) calcLMatrix();
            else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
            if (S==0) calcS(NumWFPoints, rmin, rmax);
            Lock->lock();
            calcEK(E, K, Points(rmin, rmax, NumWFPoints), el_S, IsoT, IsoF, NumFreePar, SMap, S, NAdCorr, TAdCorr, PAdCorr,
                   RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr, NumFreeAdCorr, Omega, q_e, q_f, SpinRRm, SpinRb, 
                   NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2, mJ, mv, Mv, NumFreePar, NEQ, QN, Sig, EQS, EObs, EDiff, minR, maxR, NumWFPoints,
                   MaxErr);
            Lock->unlock();
            if (NumWantedValues > 0) fillWantedValueEQS(NEQ, Sig, EQS, EDiff);
            FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints, false, 0, 0, 0, 0, (robustWeighting ? Sig : 0), 1.0, hCi_cP);
            if (done == 4 || BadCount == 20) 
            {
                if (adjustTEWeightFact)
                {
                    if (FQS < aFQS) for (m=e=0; e < nEs; e++) for (n=0; n < NEL[e]; n++) 
                        WeightF[m++] = (fabs(EL[e][n].DevR) < 4.0 ? (fabs(EL[e][n].DevR) < 1.0 ? 0.5 : 1.0) : 2.0);
                    for (n=0, done=4; n<m; n++) if (WeightF[n] == 2.0) done++;
                    if (done < 5) break;
                    for (m=e=0; e < nEs; e++) for (n=0; n < NEL[e]; n++, m++) 
                    {
                        EL[e][n].WeightFactor *= WeightF[m];
                        if (WeightF[m] != 1.0) printf("WeightFactor level with iso=%d, v=%d, J=%d now %g\n", 
                                                   EL[e][n].Iso, EL[e][n].v, EL[e][n].J, EL[e][n].WeightFactor);
                    }
                    for (n=0; n<m; n++) Sig[n] *= sqrt(WeightF[n]);
                }
                else break;
            }
            if (abs(FQS - aFQS) <  1e-10 * aFQS || (threshhold != 1e-10 && FQS < aFQS && aFQS - FQS < threshhold)) done++;
            if (FQS < aFQS && FQS > 0.0 && !isnan(FQS))
            {
                calcFitResult();
                emit potentialImproved(getPotentialData());
                if (done > 4) 
                {
                    lambda = 1e-3;
                    if (aFQS != 9e99) aFQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints);
                    else aFQS = FQS;
                    done = 0;
                }
                else 
                {
                    lambda *= 0.1;
                    aFQS = FQS;
                }
                BadCount = 0;
                calcLMEQS(NEQ + NumWantedValues, EQS, Sig, EDiff, LMEQS, hCi_cP);
                saveCoefficients(bC);
            }
            else 
            {
                restorePotential(bC);
                BadCount++;
                if (done > 4) aFQS = 9e99;
                else lambda *= 10.0;
            }
        }
        if (FQS > aFQS || FQS == 0.0 || isnan(FQS))
        {
            restorePotential(bC);
            FQS = aFQS;
        }
        Destroy(LMEQS, NumFreePar);
        Destroy(tLMEQS, NumFreePar);
    }
    else calcFitResult();
    for (F=0; F < NFC; F++)
    {
        for (I=0; I < NIso; I++)
        {
            for (J=0; J <= mJ[I]; J++) if (MaxErr[F][I][J] != 0) delete[] MaxErr[F][I][J];
            if (MaxErr[F][I] != 0) delete[] MaxErr[F][I];
        }
        delete[] MaxErr[F];
    }
    delete[] MaxErr;
    if (adjustTEWeightFact)
    {
        if (WeightFRow != 0) delete[] WeightFRow;
        WeightFRow = new int[NLevels];
        for (e=m=0, WeightSum = 0.0; e < nEs; e++) for (n=0; n < NEL[e]; n++)
        {
            WeightSum += (WeightF[m] = EL[e][n].WeightFactor);
            Row[m++] = EL[e][n].Row;
        }
        NWeightF = m;
    }
    else for (e=m=0, WeightSum = 0.0; e < nEs; e++) for (n=0; n < NEL[e]; n++) WeightSum += EL[e][n].WeightFactor;
    for (e=0; e < nGS; e++) for (n=0; n < nGL[e]; n++) WeightSum += GL[e][n].WeightFact;
    DestroyFitData(E, K, mJ, mv, NEQ, EQS, EObs, EDiff, Sig, Err, QN);
    destroyFitData();
    delete[] bC;
    return FQS;
}

void PotWorker::getSplineMinMaxR(double& Min, double& Max)
{
    if (points != 0)
    {
        Min = points[0].x;
        Max = points[numSplinePoints - 1].x;
    }
    else Min = Max = 0.0;
}

bool PotWorker::isStopFit()
{
    QMutexLocker Locker(Lock);
    if (!StopFit) return false;
    FitStopped = true;
    emit fitStopped();
    return true;
}

bool PotWorker::isFitStopped()
{
    printf("Lock.lock() PotWorker::isFitStopped\n");
    Lock->lock();
    return FitStopped;
}

void PotWorker::calc_hCi_cP_EQS(int NEQ, double** EQS, double* Sig, double* EDiff, double hCi_cP)
{
    int n, m, l, i, j;
    double CC;
    for (n=i=0; i < NLRC; i++) if (PLRC[i] >= 12 && LRCfree[i]) n++; 
    for (m=i=0; i < NLRC; i++) if (PLRC[i] >= 12 && LRCfree[i])
    {
        Sig[j = NEQ - n + m++] = 1.0 / ((CC = pow(LRC[i-1] / LRC[i-2], 3.0) * LRC[i-3]) * hCi_cP);
        EDiff[j] = CC - LRC[i];
        for (l=0; l < NumFreePar; l++) EQS[j][l] = 0.0;
        for (l=0; (l < NumFreePar ? SMap[l] - numSplinePoints - 1 < i-3 : false); l++) ;
        if (SMap[l] - numSplinePoints == i-2) EQS[j][l++] -= CC / LRC[i-3];
        if (SMap[l] - numSplinePoints == i-1) EQS[j][l++] += 3.0 * CC / LRC[i-2];
        if (SMap[l] - numSplinePoints == i) EQS[j][l++] -= 3.0 * CC / LRC[i-1];
        if (SMap[l] - numSplinePoints - 1 == i) EQS[j][l] = 1.0;
    }
}

void PotWorker::calcLMEQS(int NEQ, double** EQS, double* Sig, double* EDiff, double** LMEQS, double hCi_cP)
{
    if (hCi_cP > 0.0) calc_hCi_cP_EQS(NEQ, EQS, Sig, EDiff, hCi_cP);
    CalcLMEQS(NEQ, EQS, Sig, EDiff, LMEQS, hCi_cP, NumFreePar);
}

double PotWorker::TTFit(double****, double*****, int, int**, double*, double**, double*,
                        double*, double****, double*, double&, double)
{
    /*int n;
    double FQS = 8e99, aFQS =9e99, minR, maxR;
    if (b != 0)
    {
        TTPot->setb(b);
        Destroy(S, NumFreePar);
        S = 0;
    }
    while (FQS < 0.9999 * aFQS)
    {
        aFQS = FQS;
        calcEK(E, K, mJ, mv, Mv, NumFreePar, NEQ, QN, Sig, EQS, EObs, EDiff, minR, maxR, MaxErr);
        //printf("Nach calcEK\n");
        FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);
        if (FQSv == -1.0) FQSv = FQS;
        SvdFit(EQS, C, EDiff, Sig, NumFreePar, NEQ, 1e-7);
        TTPot->updateCoefficients(C);
        for (n=0; n < NAdCorr; n++) adCorr[n] += (C[NumFreePar - NAdCorr + n] * AdCorrWF[n]);
    }
    return FQS;*/
    return 0.0;
}

double PotWorker::FitTangToenniesPot()
{
    /*int n, NC = TTPot->getNFitCoefficients(), NLines = 0, NEQ, **QN, F, I, J;
    double ****MaxErr, ****E, *****K, **EQS, *EObs, *EDiff, *Sig, *Err;
    double FQS = 8e99, FQSv = -1.0, NFC = int(2.0 * el_S) + 1;
    double C[NumFreePar = NC + NAdCorr];
    if (SMap == 0)
    {
        SMap = new int[NumFreePar];
        for (n=0; n < NumFreePar; n++) SMap[n] = n;
    }
    estMinErr(MaxErr, NFC, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
    for (n=0; n < nEs; n++) NLines += NEL[n];
    for (n=0; n < nGS; n++) NLines += nGL[n];
    if (nEs + nGS == 0) return 0.0; 
    CreateEK(E, K, mJ, mv, NumFreePar);
    CreateFitEQS(NEQ, EQS, EObs, EDiff, Sig, Err, QN, NumFreePar, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);
    //double b1 = 0.0, b2 = TTPot->getb(), b3 = b2 + 0.1;
    FQS = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv); //FQS2*/
    
    /*double FQS1, FQS3 = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b3);
    while (FQS3 < FQS2 && b3 < 4.5)
    {
        b1 = b2;
        FQS1 = FQS2;
        b2 = b3;
        FQS2 = FQS3;
        FQS3 = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b3 += 0.1);
        }
        if (b1 == 0.0) 
        {
            FQS1 = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b1 = b2 - 0.1);
            while (FQS1 < FQS2 && b1 >= 0.4)
            {
            b3 = b2;
            FQS3 = FQS2;
            b2 = b1;
            FQS2 = FQS1;
            FQS1 = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b1 -= 0.1);
        }
    }
    ParabInterpol(b1, FQS1, b2, FQS2, b3, FQS3, b, FQS);
    for (n=0; n < 25 && fabs(FQS2 - FQS) >= 1e-8; n++)
    {
        if (b > b2)
        {
            b1 = b2;
            FQS1 = FQS2;
        }
        else 
        {
            b3 = b2;
            FQS3 = FQS2;
        }
        FQS2 = TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b2 = b);
        ParabInterpol(b1, FQS1, b2, FQS2, b3, FQS3, b, FQS);
    }
    TTFit(E, K, NEQ, QN, Sig, EQS, EObs, EDiff, MaxErr, C, FQSv, b);
    calcE(E, mJ, mv, Mv, MaxErr);
    FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);*/
    
    /*for (F=0; F < NFC; F++)
    {
        for (I=0; I < NIso; I++)
        {
            for (J=0; J <= mJ[I]; J++) if (MaxErr[F][I][J] != 0) delete[] MaxErr[F][I][J];
            if (MaxErr[F][I] != 0) delete[] MaxErr[F][I];
        }
        delete[] MaxErr[F];
    }
    delete[] MaxErr;
    DestroyFitData(E, K, mJ, mv, NEQ, EQS, EObs, EDiff, Sig, Err, QN);
    destroyFitData();*/
    return 0.0; // FQS;
}

double PotWorker::FitIsoMass(int FIso, double ISS, double &RMass, int NumWFPoints)
{
    int n, NLines = 0, J, aNAdcorr = NAdCorr, F, NFC = int(2.0 * el_S) + 1;
    double ****MaxErr, ****E, aIsoF = IsoF[FIso];
    double FQS = 8e99, FQSv = -1.0, UFact = 2e-18 * C_c * C_h * C_u / (C_hq * C_hq);
    NAdCorr = 0;
    estMinErr(MaxErr, NFC, mJ, mv, nEs, NEL, EL, nGS, nGL, GL);
    for (n=0; n < nEs; n++) NLines += NEL[n];
    for (n=0; n < nGS; n++) NLines += nGL[n];
    if (nEs + nGS == 0) return 0.0; 
    E = new double***[NFC];
    for (F=0; F < NFC; F++)
    {
        E[F] = new double**[NIso];
        E[F][FIso] = new double*[mJ[FIso] + 1];
        for (J=0; J <= mJ[FIso]; J++)
        {
            if (mv[FIso][J] >= 0) E[F][FIso][J] = new double[mv[FIso][J] + 1];
            else E[F][FIso][J] = 0;
        }
    }
    double M1 = 0.0, M2 = IsoT->redMass[FIso], M3 = M2 + ISS;
    double FQS2 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M2, NumWFPoints);
    double FQS1, FQS3 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M3, NumWFPoints);
    while (FQS3 < FQS2)
    {
        M1 = M2;
        FQS1 = FQS2;
        M2 = M3;
        FQS2 = FQS3;
        FQS3 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M3 += ISS, NumWFPoints);
    }
    if (M1 == 0.0) 
    {
        FQS1 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M1 = M2 - ISS, NumWFPoints);
        while (FQS1 < FQS2)
        {
            M3 = M2;
            FQS3 = FQS2;
            M2 = M1;
            FQS2 = FQS1;
            FQS1 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M1 -= ISS, NumWFPoints);
        }
    }
    ParabInterpol(M1, FQS1, M2, FQS2, M3, FQS3, RMass, FQS);
    for (n=0; n < 25 && fabs(FQS2 - FQS) >= 1e-8; n++)
    {
        if (RMass > M2)
        {
            M1 = M2;
            FQS1 = FQS2;
        }
        else 
        {
            M3 = M2;
            FQS3 = FQS2;
        }
        FQS2 = IsoMFit(E, MaxErr, FQSv, FIso, UFact, M2 = RMass, NumWFPoints);
        ParabInterpol(M1, FQS1, M2, FQS2, M3, FQS3, RMass, FQS);
    }
    FQS = IsoMFit(E, MaxErr, FQSv, FIso, UFact, RMass, NumWFPoints);
    for (F=0; F < NFC; F++)
    {
        for (J=0; J <= mJ[FIso]; J++) if (MaxErr[F][FIso][J] != 0)
        {
            delete[] MaxErr[F][FIso][J];
            delete[] E[F][FIso][J];
        }
        if (MaxErr[F][FIso] != 0) delete[] MaxErr[F][FIso];
        delete[] E[F][FIso];
        delete[] MaxErr[F];
        delete[] E[F];
    }
    delete[] MaxErr;
    delete[] E;
    destroyFitData();
    NAdCorr = aNAdcorr;
    IsoF[FIso] = aIsoF;
    return FQS;
}

double PotWorker::IsoMFit(double**** E, double**** MaxErr, double& FQSv, int Iso, double UFakt, double Mass, int NumWFPoints)
{
    IsoF[Iso] = UFakt * Mass;
    calcE(E, mJ, mv, Mv, NumWFPoints, MaxErr);
    double FQS = calcFQS(E, nEs, NEL, EL, nGS, nGP, nGL, PS, GL, NumWFPoints);
    if (FQSv == -1.0) FQSv = FQS;
    return FQS;
}

bool PotWorker::setMolData(IsoTab* nIsoT, double* nMU, double* nIsoF)
{
    if (Fit->isRunning()) return false;
    if (IsoT != 0) delete IsoT;
    if (NIso > 0)
    {
        delete[] MU;
        delete[] IsoF;
    }
    IsoT = nIsoT;
    NIso = IsoT->numIso;
    MU = nMU;
    IsoF = nIsoF;
    return true;
}

void PotWorker::takeMolData(IsoTab*& rIsoT, double*& rMU, double*& rIsoF)
{
    if (Fit->isRunning())
    {
        rIsoT = 0;
        rMU = 0;
        rIsoF = 0;
        return;
    }
    rIsoT = IsoT;
    IsoT = 0;
    rMU = MU;
    MU = 0;
    rIsoF = IsoF;
    IsoF = 0;
    NIso = 0;
}

void PotWorker::getSpinRGamma(int& rNSpinRGamma, double*& rSpinRGamma, double& rSpinRR1, double& rSpinRR2, double& rSpinRb, 
                              double& rSpinRRm)
{
    rNSpinRGamma = NSpinRGamma;
    rSpinRGamma = SpinRGamma;
    rSpinRR1 = SpinRR1;
    rSpinRR2 = SpinRR2;
    rSpinRb = SpinRb;
    rSpinRRm = SpinRRm;
}

void PotWorker::setSpinRGamma(int nNSpinRGamma, double* nSpinRGamma, double nSpinRR1, double nSpinRR2, double nSpinRb, double nSpinRRm)
{
    NSpinRGamma = nNSpinRGamma;
    if (SpinRGamma != 0) delete[] SpinRGamma;
    SpinRGamma = nSpinRGamma;
    SpinRR1 = nSpinRR1;
    SpinRR2 = nSpinRR2;
    SpinRb = nSpinRb;
    SpinRRm = nSpinRRm;
}

void PotWorker::cdConnectLR1C()
{
    if (Fit->isRunning()) return;
    if (points != 0)
    {
        int n;
        double dx = points[numSplinePoints - 1].x - points[numSplinePoints - 2].x, Lys = 0.0, dE = Uinf - points[numSplinePoints - 1].y;
        double ys = (points[numSplinePoints - 1].y - points[numSplinePoints - 2].y) / dx 
                  + 1.0/6.0 * dx * (points[numSplinePoints - 2].yss + 2.0 * points[numSplinePoints - 1].yss);
        for (n=0; n < NLRC - 1; n++) Lys += double(PLRC[n]) * LRC[n] * pow(points[numSplinePoints - 1].x, -PLRC[n] - 1);
        LRC[NLRC - 1] = (ys - Lys) / double(PLRC[NLRC - 1]) * pow(points[numSplinePoints - 1].x, PLRC[NLRC - 1] + 1);
        for (n=0; n < NLRC; n++) dE -= LRC[n] * pow(points[numSplinePoints - 1].x, -PLRC[n]);
        for (n=0; n < numSplinePoints; n++) points[n].y += dE;
        iO += dE;
    }
}

void PotWorker::cdConnectLR(int p1)
{
    if (Fit->isRunning()) return;
    if (numSplinePoints >= 2)
    {
        int n, i, p2;
        double yd = points[numSplinePoints - 1].y, Ra = points[numSplinePoints - 1].x, Rd = Ra - points[numSplinePoints - 2].x, B;
        double ysd = (yd - points[numSplinePoints - 2].y) / Rd + 1.0 / 6.0 * Rd * points[numSplinePoints - 2].yss 
                                                               + 1.0 / 3.0 * Rd * points[numSplinePoints - 1].yss; 
        for (i=0; (i < NLRC ? PLRC[i] != p1 : false); i++) ;
        if (i+1 < NLRC) p2 = PLRC[i+1];
        else
        {
            double *nLRC = new double[i+2];
            int *nPLRC = new int[i+2];
            bool *nLRCfree = new bool[i+2];
            p2 = p1 + 2;
            if (NLRC > 0)
            {
                for (n=0; n < NLRC; n++)
                {
                    nLRC[n] = LRC[n];
                    nPLRC[n] = PLRC[n];
                    nLRCfree[n] = LRCfree[n];
                }
                delete[] LRC;
                delete[] PLRC;
                delete[] LRCfree;
            }
            LRC = nLRC;
            PLRC = nPLRC;
            LRCfree = nLRCfree;
            NLRC = i+2;
        }
        for (n=0, yd -= Uinf; n < NLRC; n++) if (n!=i && n!=i+1)
        {
            yd += (B = LRC[n] * pow(Ra, -PLRC[n]));
            ysd -= B * double(PLRC[n]) / Ra;
        }
        LRC[i+1] = pow(Ra, p2) * (B = -(double(p1) * yd + Ra * ysd) / double(p1 - p2));
        LRC[i] = -pow(Ra, p1) * (yd + B);
        PLRC[i] = p1;
        PLRC[i+1] = p2;
        LRCfree[i] = LRCfree[i+1] = true;
    }
}

void PotWorker::cdConnectSR()
{
    if (Fit->isRunning()) return;
    if (points != 0)
    {
        double xd = points[1].x - points[0].x, ds = 0.16666666666667;
        iA = ds * pow(points[0].x, iExp + 1.0) 
            * (ds * xd * points[1].yss - (points[1].y - points[0].y) / xd);
        if (iA < 0.0)
        {
            points[0].y += points[1].y - points[2].y;
            iA = ds * pow(points[0].x, iExp + 1.0) 
                * (ds * xd * points[1].yss - (points[1].y - points[0].y) / xd);
        }
        iO = points[0].y - iA * pow(points[0].x, -iExp);
    }
}

void PotWorker::calcLMatrix()
{
    int N = numSplinePoints;
    //printf("calcLMatrix: numSplinePoints=%d\n", N);
    int st = N - 1, M;
    if (st==0)
    {
        printf("Potential::calcLMatrix: error, not enough spline points defined!\n");
        return;
    }
    if (L!=0) Destroy(L, numSplinePoints);
    L = Create(N, M=st+3);
    
    CalcLMatrix(L, points, iExp, N, Ro != 0.0);
        
    /*int i, n;  // Uncomment this block and comment the previous line if you don't want to use code from Numerical Recipes
    double B[st], F, dx[st];
    
    //for (n=0; n<N; n++) 
    //{
    //    for (i=0; i<= N+1; i++) printf("L[%d][%d]=%g ", n, i, L[n][i]);
    //    printf("\n");
    //}
    
    for (i=0; i < st; i++) dx[i] = points[i+1].x - points[i].x;
    for (n=0; n <= st; n++) for (i=0; i<M; i++) L[n][i] = 0.0;
    F = 0.5;
    B[0] = 1.0 / (points[2].x - points[0].x - F * 0.5 * dx[0]);
    L[0][0] = 3.0 * iExp * pow(points[0].x, -iExp - 1.0);
    L[0][1] = -3.0 / dx[0];
    L[0][2] = -L[0][1];
    for (n=1; n < st; n++) 
    {
        L[n][n] = L[n-1][n+1];
        L[n][n+2] = 3.0 / dx[n];
        L[n][n+1] = - L[n][n] - L[n][n+2];
    }
    L[n][n] = -1.0 * L[n-1][n+1];
    L[n][n+1] = L[n-1][n+1];
    L[n][n+2] = -3.0;
    for (i=0; i<3; i++) L[1][i] -= F * L[0][i];
    for (n=1; n < st - 1; n++)
    {
        F = 0.5 * dx[n] * B[n-1];
        B[n] = 1.0 / (points[n+2].x - points[n].x - F * 0.5 * dx[n]);
        for (i=0; i<n+3; i++) L[n+1][i] -= F * L[n][i];
    }
    F = -0.5 * dx[n] * B[n-1];
    B[n] = 1.0 / ((-1.0 - F * 0.5) * dx[n]);
    for (i=0; i<n+3; i++) L[n+1][i] -= F * L[n][i];
    for (i=0; i<M; i++) L[st][i] *= B[st -1];
        //        (L[st][i] - 0.5 * dx[st + 1] * points[st+2].yss) * B[st - 1];
    for (n = st - 2; n>=0; n--) for (i=0; i<M; i++) 
            L[n+1][i] = (L[n+1][i] - 0.5 * dx[n+1] * L[n+2][i]) * B[n];
    for (i=0; i<M; i++) L[0][i] = (L[0][i] - 0.5 * dx[0] * L[1][i]) / dx[0];
    
    //for (n=0; n<N; n++) 
    //{
    //    for (i=0; i<= N+1; i++) printf("L[%d][%d]=%g ", n, i, L[n][i]);
    //    printf("\n");
    //}
    
    //double S = 0.0;
    //for (n=0; n < NLRC; n++) S += double(PLRC[n]) * LRC[n] * pow(points[N-1].x, -1 - PLRC[n]);
    */
}

void PotWorker::calcS(int N, double RMin, double RMax)
{
    int SN = (NLRC == 0 ? 1 : 0) /*(State != 0 ? State->getStateNum() : 1)*/, NL = numSplinePoints + 2;
    int i=0, j, k, K = NL, c;
    //printf("Potential::calcS: K=%d, N=%d\n", K, N);
    /*if (SLIC6) 
    {
        K--;
        if (Uinf == 0.0 && points[numSplinePoints - 1].x > 100.0) Uinf = points[numSplinePoints - 1].y;
    }*/
    double d = (RMax - RMin) / (N-1), xd, R = RMin, R0 = points[numSplinePoints - 1].x, A, B, C, D;
    double F[NLRC], f = 1.0 / 6.0, Off[NLRC], x, iF = 0.0;
    if (RMin >= RMax) 
    {
        printf("Potential::calcS error: invalid radii!\n");
        return;
    }
    //if (NLRC == 0) SN = 1;
    if (SN == 0)
    {
        K = numSplinePoints + NLRC + 1;
        for (i=1, F[0] = pow(R0, j = -PLRC[0] - 1); i < NLRC; i++) 
            for (F[i] = F[i-1]; j >= -PLRC[i]; j--, F[i] /= R0) ;
        for (c=0; c < NLRC; c++) 
        {
            Off[c] = 1.0 / (R0 * F[c]);
            F[c] *= PLRC[c];
        }
    }
    if (NAdCorr > 0) K += NAdCorr;
    if (S!=0) Destroy(S, NumS);
    S = Create(NumS = K, N);
    SRMin = RMin;
    SRMax = RMax;
    //printf("Vor Schleife\n");
    if (R < points[0].x)
    {
        A = pow(points[0].x, -iExp);
        for (i=0; R < points[0].x; R+=d, i++) 
        {
            //printf("R=%f\n", R);
            S[0][i] = (pow(R, -iExp) - A) * (iF = 1.0 / A);
            S[1][i] = 1.0;
            for (k=2; k<K; k++) S[k][i] = 0.0;
            if (SN == 0) for (k=0, S[numSplinePoints][i] = -1.0; k < NLRC; k++) 
                S[numSplinePoints + k + 1][i] = -1.0;
        }
    }
    //printf("Vor 2.\n");
    for (j=1; i<N; i++, R+=d)
    {
        //printf("i=%d, j=%d\n", i, j);
        if (R > points[j].x) j++;
        if (j == numSplinePoints) break;
        xd = points[j].x - points[j-1].x;
        A = (points[j].x - R) / xd;
        B = 1.0 - A;
        C = f*A * (A*A - 1.0) * xd * xd;
        D = f*B * (B*B - 1.0) * xd * xd;
        for (k=0; k < NL; k++) S[k][i] = C * L[j-1][k] + D * L[j][k];
        S[j][i] += A;
        S[j+1][i] += B;
        S[0][i] *= iF;
        if (SN == 0)
        {
            for (c=1; c < NLRC; c++, k++) S[k][i] = S[NL - 1][i] * F[c] * Off[c] - 1.0;
            S[NL - 1][i] *= (F[0] * Off[0]);
            S[NL - 1][i] -= 1.0;
            S[numSplinePoints][i] -= 1.0;
        }
    }
    //printf("Nach Schleife\n");
    if (i<N)
    {
        if (SN == 0)
        {
            for ( ; i<N; R+=d, i++)
            {
                for (j=0; j < NL - 1; j++) S[j][i] = 0.0;
                //S[NL - 2][i] = 1.0;
                for (S[j++][i] = -pow(R0 = 1.0 / R, k = PLRC[0]), c=1; c < NLRC; c++, j++)
                    for (S[j][i] = S[j-1][i]; k < PLRC[c]; S[j][i] *= R0, k++) ;
                for (j = NL - 1, c=0; c < NLRC; c++, j++) S[j][i] *= Off[c]; 
            }
        }
        else
        {
            A = 4.0 * (f = pow(R0, 6.0));
            B = 0.5 * (f *= R0);
            C = -3.0 * (f *= R0);
            D = -0.5 * (f *= R0);
            for ( ; i<N; R+=d, i++) 
            {
                for (j=0; j<K-2; j++) S[j][i] = 0.0;
                f = pow(R0 = 1.0 / R, 6.0);
                S[K-2][i] = A*f;
                S[K-1][i] = B*f;
                S[K-2][i] += C * (f *= R0 * R0);
                S[K-1][i] += D*f;
            }
        }
    }
    if (NAdCorr > 0) for (R = RMin, i=0; i<N; i++, R+=d)
    {
        S[j = K - NAdCorr][i] = pow(2.0 * AdCorr_Rm / (R + AdCorr_Rm), PAdCorr);
        if (NAdCorr > 1) for (c=1, j++, x = (R - AdCorr_Rm) / (R + AdCorr_b * AdCorr_Rm); c < NAdCorr; c++, j++)
            S[j][i] = S[j-1][i] * x;
    }
    
    
    /*printf("Vor Test\n");
    double St = 0.0;
    R0 = points[numSplinePoints - 1].x;
    R = RMin + 5.0 * xd;
    d*=1000.0;
    double y;
    for (i=0, R = RMin; i<N; i+=1000, R+=d)
    {
        for (k=1, y = iA * S[0][i]; k<K-1; k++) y += S[k][i] * points[k-1].y;
        y += St * S[k][i];
        if (R > R0) 
        {
            A = R0 / R;
            y += (1.0 + pow(A, 6.0) * (3.0 * A * A - 4.0)) * Uinf;
        }
        printf("R=%f, y=%f=%f\n", R, getPoint(R), y);
    }*/
}

void PotWorker::getS(double &Rmin, double &Rmax, int &NPoints, int& NC, double**& rS)
{
    if (Type == SplinePotential && L==0) calcLMatrix();
    else if (Type != SplinePotential && NAdCorr > 0 && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
    calcS(NPoints, Rmin, Rmax);
    NC = NumFreePar;
    rS = S;
    S=0;
}

void PotWorker::DestroyLS()
{
    if (L!=0) 
    {
        Destroy(L, numSplinePoints);
        L=0;
    }
    if (S!=0)
    {
        Destroy(S, NumS);
        S=0;
    }
    if (SMap != 0)
    {
        delete[] SMap;
        SMap = 0;
    }
}

double PotWorker::getMinR(double E, double Prec)
{
    Lock->lock();
    double R, R1 = rmin, R2 = R1 + 0.1, E1 = Point(R1), E2 = Point(R2);
    int n;
    for (n=0; abs(E2 - E) > Prec && n < 1000; n++)
    {
        R = R2 + (E2 - E) / (E1 - E2) * (R2 - R1);
        R1 = R2;
        E1 = E2;
        R2 = R;
        E2 = Point(R);
    }
    Lock->unlock();
    printf("getMinR: R=%f, n=%d\n", R2, n);
    if (n==1000) return -1.0;
    return R2;
}

double *PotWorker::getPoints(double Ri, double Ra, int N, int FC)
{
    //printf("Potential::getPoints, Ri=%f, Ra=%f, N=%d\n", Ri, Ra, N);
    QMutexLocker lock(Lock);
    return Points(Ri, Ra, N, FC);
}

double *PotWorker::get_dVdR(const double Rmin, const double Rmax, const int numPoints) const
{
    QMutexLocker lock(Lock);
    return dVdR(Rmin, Rmax, numPoints);
}

PotentialData* PotWorker::getPotentialData()
{
    int n;
    PotentialData *R = new PotentialData();
    Lock->lock();
    R->N = numSplinePoints;
    R->NLRC = NLRC;
    if (NLRC > 0)
    {
        R->PLRC = new int[NLRC];
        R->LRC = new double[NLRC];
        R->LRCFree = new bool[NLRC];
        for (n=0; n < NLRC; n++)
        {
            R->PLRC[n] = PLRC[n];
            R->LRC[n] = LRC[n];
            R->LRCFree[n] = LRCfree[n];
        }
    }
    R->Ro = Ro;
    R->NAdCorr = NAdCorr;
    R->TAdCorr = TAdCorr;
    R->PAdCorr = PAdCorr;
    R->points = new SplinePoint[numSplinePoints];
    for (n=0; n < numSplinePoints; n++) R->points[n] = points[n];
    R->iA = iA;
    R->iO = iO;
    R->iExp = iExp;
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
    R->AdCorr_b = AdCorr_b;
    R->AdCorr_Rm = AdCorr_Rm;
    R->Uinf = Uinf;
    R->shortRVariable = shortRVariable;
    R->fitResult = Result;
    R->FQS = FQS;
    Lock->unlock();
    return R;
}

double* PotWorker::Points(double Ri, double Ra, int N, int FC)
{
    if (Type == NoPotential) return 0;
    double *R = new double[N], r, s = (Ra - Ri) / double(N - 1);
    int n;
    for (n=0, r = Ri; n<N; n++, r+=s) R[n] = Point(r, FC);
    return R;
}

double PotWorker::getPoint(double R, int FC)
{
    Lock->lock();
    double U = Point(R, FC);
    Lock->unlock();
    return U;
}

double PotWorker::Point(double R, int /*FC*/)
{
    double U = 0.0;
    if (points == 0) return U;
    int n, N = numSplinePoints;
    for (n=0; (n < N ? points[n].x < R : false); n++) ;
    if (n==0) 
    {
        if (R < points[0].x && Ro == 0.0) 
        {
            if (iA == 0.0 || reCalcSC)
            {
                printf("getPoints: connect SR\n");
                double xd = points[1].x - points[0].x, ds = 0.16666666666667;
                iA = ds * pow(points[0].x, 7.0) 
                    * (ds * xd * points[1].yss - (points[1].y - points[0].y) / xd);
                iO = points[0].y - iA * pow(points[0].x, -6.0);
                emit potentialChanged();
                reCalcSC = false;
            }
            return U + iO + iA * pow(R, -iExp);
        }
        return U + points[0].y;
    }
    if (n==N || (Ro > 0.0 && R > Ro)) 
    {
        //printf("Long Range, R=%f, N=%d\n", R, N);
        double Ra = points[N-1].x, Rb, Res;
        if (NLRC == 0 || reCalcLC)
        {
            printf("getPoints: connect LR\n");
            if (NLRC >= 2) cdConnectLR(PLRC[NLRC - 2]);
            else
            {
                if (NLRC == 0)
                {
                    LRC = new double[NLRC = 2];
                    PLRC = new int[NLRC];
                    PLRC[0] = 6;
                    PLRC[1] = 8;
                }
                Rb = pow(Ra, PLRC[0]);
                double ys = (points[N-1].y - points[N-2].y) / (points[N-1].x - points[N-2].x)
                        + (points[N-1].x - points[N-2].x) * (points[N-2].yss + 2.0 * points[N-1].yss) / 6.0;
                LRC[0] = Rb * (0.5 * double(2 + PLRC[0]) * (Uinf - points[N-1].y) - 0.5 * (Ra) * ys);
                Rb *= Ra * Ra;
                LRC[1] = Rb * (0.5 * double(PLRC[0]) * (points[N-1].y - Uinf) + 0.5 * Ra * ys);
                emit potentialChanged();
            }
            reCalcLC = false;
        }
        //printf("y=%g, g=%g\n", points[N-1].y, Uinf - LRC[0] * pow(points[N-1].x, -1 * PLRC[0])
            //        - LRC[1] * pow(points[N-1].x, -1 * PLRC[1]));
        int p = PLRC[NLRC-1];
        for (n = NLRC - 2, Res = LRC[NLRC - 1]; n>=0 ; n--)
        {
            //printf("n=%d, Res=%g\n", n, Res);
            while (p > PLRC[n])
            {
                p--;
                Res /= R;
                //printf("R=%g, Res=%g, p=%d\n", R, Res, p);
            }
            Res += LRC[n];
            //printf("Res=%g, PLRC[%d]=%d\n", Res, n, PLRC[n]);
        }
        Res *= pow(R, -1 * PLRC[0]);
        //printf("R=%g, points[N-1].x=%g, Result=%g\n", R, points[N-1].x, Uinf - Res);
        return U + Uinf - Res;
    }
    /*if (n>=N-1 && SLIC6) 
    {
        //printf("SLIC6!\n");
        n = N-1;
        return U + points[n].y - points[n].yss * (pow(points[n].x, -6.0) - pow(R, -6.0));
    }*/
    double xd = points[n].x - points[n-1].x;
    double A = (points[n].x - R) / xd;
    double B = 1.0 - A;
    //printf("R=%f, n=%d, N=%d, xd=%f, A=%f, B=%f\n", R, n, N, xd, A, B);
    return U + A * points[n-1].y + B * points[n].y + xd * xd * (A*(A*A-1.0) * points[n-1].yss
                       + B*(B*B-1.0) * points[n].yss) * 0.16666666666667;
}

double* PotWorker::dVdR(const double Rmin, const double Rmax, const int numPoints) const
{
    const double h = (Rmin - Rmax) / (numPoints - 1);
    double *Ret = new double[numPoints], R = rmin;
    int n = 0;
    if (Rmin < points[0].x)
    {
        const double F = -iExp * iA, exp = -iExp - 1.0;
        for ( ; R < points[0].x && n < numPoints; ++n, R+=h) Ret[n] = F * pow(R, exp);
    }
    if (Rmin < points[numSplinePoints - 1].x && n < numPoints)
    {
        int p = 1;
        double deltaX, dDeltaX, A, B, Q, deltaXdsix, yssF1, yssF2;
        const double one = 1.0, three = 3.0, dsix = one / 6.0;
        for ( ; R < points[numSplinePoints - 1].x && n < numPoints; ++n, R+=h)
        {
            if (R > points[p].x)
            {
                while (R > points[p].x) ++p;
                dDeltaX = 1.0 / (deltaX = points[p].x - points[p-1].x);
                Q = (points[p].y - points[p-1].y) * dDeltaX;
                deltaXdsix = deltaX * dsix;
                yssF1 = deltaXdsix * points[p-1].yss;
                yssF2 = deltaXdsix * points[p].yss;
            }
            B = 1.0 - (A = (points[p].x - R) * dDeltaX);
            Ret[n] = Q - (three * A * A - one) * yssF1 + (three * B * B - one) * yssF2;
        }
    }
    if (n < numPoints)
    {
        double Res, dR;
        const double one = 1.0;
        int p = PLRC[NLRC-1], i;
        for (i = NLRC - 2, Res = LRC[NLRC - 1]; i>=0 ; --i)
        {
            dR = one / R;
            while (p > PLRC[i])
            {
                --p;
                Res *= dR;
            }
            Res += static_cast<double>(PLRC[i]) * LRC[i];
        }
        Ret[n] = Res * pow(dR, one + static_cast<double>(PLRC[0]));
    }
    return Ret;
}

void PotWorker::getMaximum(double &R, double &U)
{
    Lock->lock();
    if (points != 0)
    {
        int n, N = numSplinePoints;
        double a, b, c, xd, yssd, r, E, x;
        for (n=0, U = points[0].y, R = points[0].x; n<N-1; n++)
        {
            yssd = points[n+1].yss - points[n].yss;
            a = (points[n+1].x * points[n].yss - points[n].x * points[n+1].yss) / yssd;
            xd = points[n+1].x - points[n].x;
            b = (2.0 * (points[n+1].y - points[n].y) + points[n].x * points[n].x * points[n+1].yss
                     - points[n+1].x * points[n+1].x * points[n].yss) / yssd - xd * xd / 3;
            if ((c=a*a-b) <= 0.0) continue;
            x = (r = sqrt(c)) - a;
            if (x >= points[n].x && x < points[n+1].x && (E = Point(x)) > U)
            {
                U = E;
                R = x;
            }
            x = - r - a;
            if (x >= points[n].x && x < points[n+1].x && (E = Point(x)) > U)
            {
                U = E;
                R = x;
            }
        }
        //printf("Potential::getMinimum: R=%f, U=%f\n", R, U);
    }
    else U = Point(R = rmin);
    Lock->unlock();
}

void PotWorker::getMinimum(double &R, double &U)
{
    Lock->lock();
    if (points != 0)
    {
        int n, N = numSplinePoints;
        double a, b, c, xd, yssd, r, E, x;
        for (n=0, U = points[0].y; n<N-1; n++)
        {
            yssd = points[n+1].yss - points[n].yss;
            a = (points[n+1].x * points[n].yss - points[n].x * points[n+1].yss) / yssd;
            xd = points[n+1].x - points[n].x;
            b = (2.0 * (points[n+1].y - points[n].y) + points[n].x * points[n].x * points[n+1].yss
                     - points[n+1].x * points[n+1].x * points[n].yss) / yssd - xd * xd / 3;
            if ((c=a*a-b) <= 0.0) continue;
            x = (r = sqrt(c)) - a;
            if (x >= points[n].x && x < points[n+1].x && (E = Point(x)) < U)
            {
                U = E;
                R = x;
            }
            //printf("p[%d].x=%f, p[%d].y=%f, x1=%f, E1=%f, ", n, points[n].x, n, points[n].y, x, E);
            x = - r - a;
            if (x >= points[n].x && x < points[n+1].x && (E = Point(x)) < U)
            {
                U = E;
                R = x;
            }
            //printf("x2=%f, E2=%f\n", x, E);
        }
        //printf("Potential::getMinimum: R=%f, U=%f\n", R, U);
    }
    /*else if (points != 0 && aPot != 0)
    {
        double x1, y1, x2, y2=0.0, x3, y3=Point(rmin), r;
        for (r = rmin + 0.1, U = y3; r <= rmax; r+=0.1)
        {
            x1 = x2;
            y1 = y2;
            x2 = x3;
            y2 = y3;
            x3 = r;
            y3 = Point(r);
            if (y1 != 0.0 && y1 >= y2 && y2 < y3 && y2 < U)
            {
                R = x2;
                U = y2;
            }
        }
        x2 = R;
        y2 = U;
        x1 = R - 0.1;
        y1 = Point(x1);
        x3 = R + 0.1;
        y3 = Point(x3);
        while (fabs(r - U) > 0.0001)
        {
            ParabInterpol(x1, y1, x2, y2, x3, y3, R, U);
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
        U = r;
    }*/
    else R = U = 0.0;
    Lock->unlock();
}

void PotWorker::shift(double E)
{
    Lock->lock();
    int n, N = numSplinePoints;
    for (n=0; n<N; n++) points[n].y += E;
    Uinf += E;
    iO += E;
    Lock->unlock();
}

bool PotWorker::addPoint(double x, double y)
{
    if (Fit->isRunning()) return false;
    if (points == 0) 
    {
        printf("Potential::addPoint error: points=0!\n");
        return false;
    }
    if (L!=0) 
    {
        Destroy(L, numSplinePoints);
        L = 0;
    }
    if (S!=0) 
    {
        Destroy(S, NumS);
        S = 0;
    }
    if (SMap != 0) 
    {
        delete[] SMap;
        SMap = 0;
    }
    int n;
    //if (aPot != 0) y -= aPot->GetPoint(x);
    SplinePoint *buff = new SplinePoint[++numSplinePoints];
    for (n=0; (n < numSplinePoints - 1 ? points[n].x < x : false); n++) buff[n] = points[n];
    buff[n].x = x;
    buff[n].y = y;
    buff[n].variable = true;
    for (n++; n < numSplinePoints; n++) buff[n] = points[n-1];
    delete[] points;
    points = buff;
    return true;
}

bool PotWorker::movePoint(int& n, double x, double y)
{
    if (Fit->isRunning()) return false;
    if (points == 0)
    {
        printf("Potential::movePoint error: points = 0!\n");
        return false;
    }
    if (n < 0 || n >= numSplinePoints)
    {
        printf("Potential::movePoint error: n=%d, numSplinePoints=%d\n", n, numSplinePoints);
        return false;
    }
    if (L!=0) 
    {
        Destroy(L, numSplinePoints);
        L = 0;
    }
    if (S!=0) 
    {
        Destroy(S, NumS);
        S = 0;
    }
    if (SMap != 0) 
    {
        delete[] SMap;
        SMap = 0;
    }
    while (n > 0 ? x < points[n-1].x : false) 
    {
        points[n] = points[n-1];
        n--;
    }
    while (n < numSplinePoints - 1 ? x > points[n+1].x : false) 
    {
        points[n] = points[n+1];
        n++;
    }
    points[n].x = x;
    points[n].y = y;
    return true;
}

bool PotWorker::removePoint(int n)
{
    if (Fit->isRunning()) return false;
    if (points == 0)
    {
        printf("Potential::removePoint error: points = 0!\n");
        return false;
    }
    if (n < 0 || n >= numSplinePoints)
    {
        printf("Potential::removePoint error: n=%d, numSplinePoints=%d\n", n, numSplinePoints);
        return false;
    }
    if (L!=0) 
    {
        Destroy(L, numSplinePoints);
        L = 0;
    }
    if (S!=0) 
    {
        Destroy(S, NumS);
        S = 0;
    }
    if (SMap != 0) 
    {
        delete[] SMap;
        SMap = 0;
    }
    for (n++; n < numSplinePoints; n++) points[n-1] = points[n];
    numSplinePoints--;
    return true;
}

SplinePoint *PotWorker::getSplinePoints(int &n)
{
    QMutexLocker lock(Lock);
    if (points == 0)
    {
        n = 0;
        return 0;
    }
    int i;
    SplinePoint *rPoints = new SplinePoint[n = numSplinePoints];
    for (i=0; i<n; i++) rPoints[i] = points[i];
    //if (aPot != 0) for (i=0; i<n; i++) rPoints[i].y -= aPot->GetPoint(points[i].x);
    return rPoints;
}

void PotWorker::calcFitResult()
{
    int s, n, N;
    double devq;
    BLLock->lock();
    for (s = N = Result.NBad = Result.NBadPAL = 0, Result.FQS_Bad = Result.Sigma = 0.0; s < nEs; s++) 
        for (n=0, N += NEL[s]; n < NEL[s]; n++)
    {
        Result.Sigma += (devq = EL[s][n].DevR * EL[s][n].DevR);
        if (fabs(EL[s][n].DevR >= 4.0))
        {
            Result.NBad++;
            Result.FQS_Bad += devq;
            if (EL[s][n].err < 1e-4) Result.NBadPAL++;
        }
    }
    for (s=0; s < nGS; s++) for (n=0, N += nGL[s]; n < nGL[s]; n++)
    {
        Result.Sigma += (devq = GL[s][n].DevR * GL[s][n].DevR);
        if (fabs(GL[s][n].DevR >= 4.0))
        {
            Result.NBad++;
            Result.FQS_Bad += devq;
            if (GL[s][n].err < 1e-4) Result.NBadPAL++;    
        }
    }
    BLLock->unlock();
    Result.Sigma = sqrt(Result.Sigma / (N - NumFreePar));
    if (NumWantedValues > 0)
    {
        for (n=1, Result.LRelParDev = fabs(WantedValues[0].Value - WantedValues[0].CurValue) / WantedValues[0].Precision; n < NumWantedValues; n++)
            if ((devq = fabs(WantedValues[n].Value - WantedValues[n].CurValue) / WantedValues[n].Precision) < Result.LRelParDev)
                Result.LRelParDev = devq;
    }
    else Result.LRelParDev = 0.0;
}
