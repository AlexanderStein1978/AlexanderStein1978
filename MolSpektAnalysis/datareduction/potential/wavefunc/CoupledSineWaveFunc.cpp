//
// C++ Implementation: CoupledSineWaveFunc
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "CoupledSineWaveFunc.h"
#include "naturalspline.h"
#include "utils.h"
#include "molecule.h"
#include "SplinePoint.h"

#include <cmath>

#include <QDialog>
#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QFile>
#include <QTextStream>


CoupledSineWaveFunc::CoupledSineWaveFunc(MainWindow *MW, Molecule *mol) 
                    : MDIChild(SimpleDiagWindow, MW, "SineWaveFuncs (*.SWF)", ".SWF")
{
    Data = 0;
    Iso = SC = 0;
    R = 0;
    nChan = nCoeff = nI = nJ = nv = 0;
    States = 0;
    E = 0;
    MC = 0;
    Mol = mol;
    ASC = 0;
    Basis = 0;
    RMin = RMax = 0.0;
    nPoints = 0;
    vA = CA = 0;
    MapFuncs = 0;
}

CoupledSineWaveFunc::~CoupledSineWaveFunc()
{
    DestroyData();
    if (MapFuncs != 0) delete MapFuncs;
}

void CoupledSineWaveFunc::DestroyData()
{
    int I, J, v, C;
    if (Data != 0)
    {
        for (I=0; I < nI; I++)
        {
            for (J=0; J < nJ; J++)
            {
                for (v=0; v < nv; v++) if (Data[I][J][v] != 0)
                {
                    for (C=0; C < nChan; C++)  delete[] Data[I][J][v][C];
                    delete[] Data[I][J][v];
                    delete[] MC[I][J][v];
                }
                delete[] Data[I][J];
                delete[] MC[I][J];
                delete[] E[I][J];
            }
            delete[] Data[I];
            delete[] MC[I];
            delete[] E[I];
        }
        delete[] Data;
        Data = 0;
        delete[] MC;
        delete[] E;
    }
    if (R != 0)
    {
        delete[] R;
        R=0;
    }
    if (Iso != 0)
    {
        delete[] Iso;
        Iso = 0;
    }
    if (SC != 0)
    {
        delete[] SC;
        delete[] States;
        SC = 0;
    }
    if (Basis != 0) 
    {
        Destroy(Basis, nPoints);
        Basis = 0;
    }
    if (vA != 0) 
    {
        Destroy(vA, nI, nJ);
        vA = 0;
    }
    if (CA != 0) 
    {
        Destroy(CA, nI, nJ);
        CA = 0;
    }
}

double**** CoupledSineWaveFunc::getEnergy()
{
    double ****rE = Create(1, nI, nv, nJ);
    int I, v, J;
    for (I=0; I < nI; I++) for (v=0; v < nv; v++) for (J=0; J < nJ; J++) rE[0][I][v][J] = E[I][J][v];
    return rE;
}

void CoupledSineWaveFunc::CalcvCA()
{
    int v, I, J, n, c[nChan], bc;
    double MMC;
    if (vA != 0) Destroy(vA, nI, nJ);
    if (CA != 0) Destroy(CA, nI, nJ);
    vA = CreateInt(nI, nJ, nv);
    CA = CreateInt(nI, nJ, nv);
    for (I=0; I < nI; I++) for (J=0; J < nJ; J++)
    {
        for (n=0; n < nChan; n++) c[n] = 0;
        for (v=0; v < nv; v++) if (Data[I][J][v] != 0)
        {
            for (n=1, MMC = MC[I][J][v][0], bc=0; n < nChan; n++) 
                if (MC[I][J][v][n] > MMC)
            {
                MMC = MC[I][J][v][n];
                bc = n;
            }
            vA[I][J][v] = c[bc]++;
            CA[I][J][v] = bc;
        }
    }
}

void CoupledSineWaveFunc::getIso(int& N, int*&I)
{
    N = nI;
    I = Iso;
}

int CoupledSineWaveFunc::getNChannels()
{
    return nChan;
}

int CoupledSineWaveFunc::getNLevels(int I, int J)
{
    if (I >= nI && J >= nJ) return 0;
    int n, r;
    for (n=r=0; r < nv; n++) if (Data[I][J][n] != 0) r++;
    return r;
}

void CoupledSineWaveFunc::getStateComp(int Channel, ElState*& State, int& Comp)
{
    State = States[Channel];
    Comp = SC[Channel];
}

void CoupledSineWaveFunc::getv(int I, int J, int& N, int*& v)
{
    if (I >= nI && J >= nJ)
    {
        N=0;
        v=0;
        return;
    }
    int n;
    for (n=N=0; n < nv; n++) if (Data[I][J][n] != 0) N++;
    v = new int[N];
    for (n=N=0; n < nv; n++) if (Data[I][J][n] != 0) v[N++] = n;
}

int CoupledSineWaveFunc::getv(int I, int J, double SE)
{
    int n;
    if (I<0 || I >= nI || J<0 || J >= nJ) return -1;
    for (n=0; n < nv; n++) if (fabs(SE - E[I][J][n]) < 1e-3) return n;
    return -1;
}

bool CoupledSineWaveFunc::isJA(int I, int J)
{
    int n;
    if (I<0 || I >= nI || J<0 || J >= nJ) return false;
    for (n=0; n < nv; n++) if (Data[I][J][n] != 0) return true;
    return false;
}

bool CoupledSineWaveFunc::getWaveFunc(int I, int J, int v, int Channel, double RStart, double RStop, int NPoints, double* WF, 
                                      double& rE, int& MMC)
{
    if (I >= nI || J >= nJ || v >= nv) return false;
    if (Data[I][J][v] == 0) return false;
    int n, m, j;
    double M;
    double h = (RStop - RStart) / (NPoints - 1), r, x, F = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1;
    //if (RStart != RMin || RStop != RMax || NPoints != nPoints) CalcBasis(RStart, RStop, NPoints); 
    for (n=0; n < NPoints; n++) WF[n] = 0.0;
    for (n=0, r = RStart; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++) 
            WF[n] += (sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                   - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)))) * Data[I][J][v][Channel][j];
        WF[n] *= F * sqrt(MapFuncs->getDerivative(r));
    }
    rE = E[I][J][v];
    for (MMC = 0, n=1, M = MC[I][J][v][0]; n < nChan; n++) if (MC[I][J][v][n] > M)
    {
        MMC = n;
        M = MC[I][J][v][n];
    }
    return true;
}

void CoupledSineWaveFunc::CalcBasis(double RStart, double RStop, int NPoints)
{
    double h = (RStop - RStart) / (NPoints - 1), r, x, F = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1;
    int n, m, j;
    if (Basis != 0) Destroy(Basis, nPoints);
    Basis = Create(nPoints = NPoints, nCoeff);
    for (n=0; n < NPoints; n++) for (j=0; j < nCoeff; j++) Basis[n][j] = 0.0;
    for (n=0, r = RMin = RStart, RMax = RStop; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++) 
            Basis[n][j] = (sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                            - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)))) * F * sqrt(MapFuncs->getDerivative(r));
    }
}

bool CoupledSineWaveFunc::getWaveFuncs(int I, int v, int *Channels, int NChannels, int NJ, int *J, 
                                       double RStart, double RStop, int NPoints, 
                                       double*** WF, double* rE, int *va, int *ca)
{
    if (I >= nI || v >= nv || J[NJ-1] > nJ) return false;
    //if (RStart != RMin || RStop != RMax || NPoints != nPoints) CalcBasis(RStart, RStop, NPoints);
    int n, m, j, i, c;
    double M;
    double h = (RStop - RStart) / (NPoints - 1), r, x, F = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1, TF;
    if (vA == 0 || CA == 0) CalcvCA();
    for (n=0; n < NPoints; n++) for (i=0; i < NJ; i++) for (c=0; c < NChannels; c++) WF[i][c][n] = 0.0;
    for (n=0, r = RStart; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++)
        {
            M = sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                            - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)));
            for (i=0; i < NJ; i++) for (c=0; c < NChannels; c++) 
                if (Data[I][J[i]][v] != 0) WF[i][c][n] += Data[I][J[i]][v][Channels[c]][j] * M;
        }
        TF = F * sqrt(MapFuncs->getDerivative(r));
        for (i=0; i < NJ; i++) for (c=0; c < NChannels; c++) WF[i][c][n] *= TF;
    }
    for (i=0; i < NJ; i++) 
    {
        rE[i] = E[I][J[i]][v];
        va[i] = vA[I][J[i]][v];
        ca[i] = CA[I][J[i]][v];
    }
    return true;
}

bool CoupledSineWaveFunc::getWaveFuncs(const int I, const int J, int &Nv, const int NChannels, const double RStart, const double RStop,
                                       const int NPoints, double ***&WFs, int*& ri, int*& ra, int*& P, double *&F, double *&N)
{
    if (I >= nI || J >= nJ) return false;
    if (Nv > nv) Nv = nv;
    int n, m, j, i, c;
    double M;
    double h = (RStop - RStart) / (NPoints - 1), r, x, Fac = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1, TF;
    WFs = Create(NChannels, Nv, NPoints);
    ri = new int[Nv];
    ra = new int[Nv];
    P = new int[Nv];
    F = new double[Nv];
    N = new double[Nv];
    if (vA == 0 || CA == 0) CalcvCA();
    for (i=0; i < Nv; i++)
    {
        for (n=0; n < NPoints; n++) for (c=0; c < NChannels; c++) WFs[c][i][n] = 0.0;
        N[i] = 0.0;
        F[i] = 1.0;
        ri[i] = P[i] = 0;
        ra[i] = NPoints - 1;
    }
    for (n=0, r = RStart; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++)
        {
            M = sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                            - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)));
            for (i=0; i < Nv; i++) for (c=0; c < NChannels; c++)
                if (Data[I][J][i] != 0) WFs[c][i][n] += Data[I][J][i][c][j] * M;
        }
        TF = Fac * sqrt(MapFuncs->getDerivative(r));
        for (i=0; i < Nv; i++) for (c=0; c < NChannels; c++)
        {
            WFs[c][i][n] *= TF;
            N[i] += WFs[c][i][n] * WFs[c][i][n];
        }
    }
    for (i=0; i < Nv; ++i) N[i] = 1.0 / N[i];
    return true;
}

bool CoupledSineWaveFunc::getWaveFuncs(int I, int J, int Channel, double RStart, double RStop, int NPoints, double** WF, 
                                       double* rE)
{
    if (I >= nI || J >= nJ) return false;
    int n, j, m, Nv, v, c;
    double M, MS;
    double h = (RStop - RStart) / (NPoints - 1), r, x, F = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1;
    //if (RStart != RMin || RStop != RMax || NPoints != nPoints) CalcBasis(RStart, RStop, NPoints);
    for (v = Nv = 0; v < nv; v++) if (Data[I][J][v] != 0) Nv = v + 1;
    if (Nv == 0) return false;
    double Sq[Nv], **WFt = Create(Nv, nChan);
    for (v=0; v < Nv; v++) Sq[v] = 0.0;
    for (n=0, r = RStart; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        for (v=0; v < Nv; v++) for (c=0; c < nChan; c++) WFt[v][c] = 0.0;
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++) for (c=0; c < nChan; c++)
        {
            M = sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                            - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)));
            for (v=0; v < Nv; v++) if (Data[I][J][v] != 0) 
                WFt[v][c] += Data[I][J][v][c][j] * M;
        }
        MS = sqrt(MapFuncs->getDerivative(r));
        for (v=0; v < Nv; v++) 
        {
            for (c=0; c < nChan; c++) 
            {
                WFt[v][c] *= F * MS;
                Sq[v] += WFt[v][c] * WFt[v][c];
            }
            WF[v][n] = WFt[v][Channel];
        }
    }
    for (v=0; v < Nv; v++) 
    {
        rE[v] = E[I][J][v];
        for (Sq[v] = 1.0 / sqrt(Sq[v]), n=0; n < NPoints; n++) WF[v][n] *= Sq[v];
    }
    Destroy(WFt, Nv);
    return true;
}

bool CoupledSineWaveFunc::getWaveFuncs(int I, int J, double RStart, double RStop, int NPoints, double*** WF, double* rE)
{
    if (I >= nI || J >= nJ) return false;
    //if (RStart != RMin || RStop != RMax || NPoints != nPoints) CalcBasis(RStart, RStop, NPoints);
    int n, m, j, Nv, v, Channel;
    double M;
    double h = (RStop - RStart) / (NPoints - 1), r, x, F = 0.5 / double(nCoeff), F1 = 0.5 * M_PI / double(nCoeff);
    double F2 = double(2 * nCoeff - 1) * F1;
    for (v = Nv = 0; v < nv; v++) if (Data[I][J][v] != 0) Nv = v + 1;
    if (Nv == 0) return false;
    for (n=0; n < NPoints; n++) for (v=0; v < Nv; v++) for (Channel = 0; Channel < nChan; Channel++)
        WF[Channel][v][n] = 0.0;
    for (n=0, r = RStart; r < R[0] && n < NPoints; r+=h, n++) ;
    for (m=0; r <= R[nCoeff - 1] && n < NPoints; r+=h, n++)
    {
        while (R[m] < r) m++;
        x = MapFuncs->getPoint(r);
        for (j=0; j < nCoeff; j++) 
        {
            M = sin(F2 * (x - double(j))) / sin(F1 * (x - double(j)))
                            - sin(F2 * (x + double(j))) / sin(F1 * (x + double(j)));
            for (v=0; v < Nv; v++) if (Data[I][J][v] != 0) for (Channel = 0; Channel < nChan; Channel++) 
                WF[Channel][v][n] += Data[I][J][v][Channel][j] * M;
        }
        for (Channel = 0; Channel < nChan; Channel++) for (v=0; v < Nv; v++) 
            WF[Channel][v][n] *= F * sqrt(MapFuncs->getDerivative(r));
    }
    for (v=0; v < Nv; v++) rE[v] = E[I][J][v];
    return true;
}

bool CoupledSineWaveFunc::readData(QString FN)
{
    if (Mol == 0 && MW != 0)
    {
        QDialog *D = new QDialog(this);
        QGridLayout *L = new QGridLayout(D);
        D->setWindowTitle("MolSpektAnalysis");
        L->addWidget(new QLabel("Please select the molecule the current potential belongs to:", D), 0, 0, 1, 2);
        L->addWidget(new QLabel("Molecule:", D), 1, 0);
        QComboBox *B = new QComboBox(D);
        int n, N = MW->getNumMolecules();
        for (n=0; n<N; ++n) B->addItem(MW->getMolecule(n)->getName());
        B->setEditable(false);
        L->addWidget(B, 1, 1);
        L->setRowMinimumHeight(2, 20);
        QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
        L->addWidget(OK, 3, 0);
        L->addWidget(Cancel, 3, 1);
        connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
        connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
        if (D->exec() == QDialog::Accepted) Mol = MW->getMolecule(B->currentIndex());
        delete D;
    }
    if (Mol == 0) 
    {
        printf("CoupledSineWaveFunc::readData error: Can't read any data without being assinged to a molecule!");
        return false;
    }
    int n, m, N, I, J, v, c;
    SplinePoint *points;
    QFile Datei(FN);
    if (!read(&Datei)) return false;
    QTextStream S(&Datei);
    QStringList L = S.readLine().split(' ', Qt::SkipEmptyParts);
    if (L.count() > 0 ? L[0].indexOf("name:", 0, Qt::CaseInsensitive) == -1 : true) return false;
    if (L.count() >= 2) setName(L[1]);
    if ((L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() > 0 ?
        L[0].indexOf("source:", 0, Qt::CaseInsensitive) == -1 : true) return false;
    if (L.count() >= 2) Source = L[1];
    DestroyData();
    if ((nChan = (L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() / 2) == 0) return false;
    for (n=0, SC = new int[nChan], States = new ElState*[nChan]; n < nChan; n++)
    {
        States[n] = Mol->getState(L[2*n]);
        SC[n] = L[2*n+1].toInt();
    }
    if ((nI = (L = S.readLine().split(' ', Qt::SkipEmptyParts)).count()) == 0) return false;
    for (n=0, Iso = new int[nI]; n < nI; n++) Iso[n] = L[n].toInt();
    if ((nCoeff = (L = S.readLine().split(' ', Qt::SkipEmptyParts)).count()) == 0) return false;
    for (n=0, R = new double[nCoeff], points = new SplinePoint[nCoeff]; n < nCoeff; n++) 
    {
        points[n].y = double(n);
        points[n].x = R[n] = L[n].toDouble();
    }
    if (MapFuncs == 0) MapFuncs = new NaturalSpline(points, nCoeff);
    else MapFuncs->setPoints(points, nCoeff);
    if ((L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() < 3) return false;
    if ((nJ = L[0].toInt()) < 0 || (nv = L[1].toInt()) < 0) return false;
    for (I=0, Data = new double****[nI], MC = new double***[nI], E = new double**[nI]; I < nI; I++)
        for (J=0, Data[I] = new double***[nJ], MC[I] = new double**[nJ], E[I] = new double*[nJ]; J < nJ; J++)
            for (v=0, Data[I][J] = new double**[nv], MC[I][J] = new double*[nv], E[I][J] = new double[nv]; v < nv; v++)
                Data[I][J][v] = 0;
    for (n=0, N = L[2].toInt(); n<N; n++)
    {
        if (S.atEnd()) break;
        if ((L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() < 4) break;
        if ((I = L[0].toInt()) >= nI || (J = L[1].toInt()) >= nJ || (v = L[2].toInt()) >= nv) break;
        if (I < 0 || J < 0 || v < 0) break;
        E[I][J][v] = L[3].toDouble();
        if ((L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() < nChan) break;
        for (c=0, MC[I][J][v] = new double[nChan], Data[I][J][v] = new double*[nChan]; c < nChan; c++)
        {
            MC[I][J][v][c] = L[c].toDouble();
            Data[I][J][v][c] = new double[nCoeff];
        }
        for (m=0; m < nCoeff; m++)
        {
            if ((L = S.readLine().split(' ', Qt::SkipEmptyParts)).count() < nChan) break;
            for (c=0; c < nChan; c++) Data[I][J][v][c][m] = L[c].toDouble();
        }
        if (m < nCoeff) break;
    }
    if (n<N) return false;
    Saved();
    return true;
}

void CoupledSineWaveFunc::setData(int NI, int NJ, int Nv, int NChan, int NCoeff, ElState **St, int *C, int* I, double* RS, 
                                  double***** data, double ****MCo, double ***En)
{
    DestroyData();
    int n;
    SplinePoint *points = new SplinePoint[NCoeff];
    for (n=0; n < NCoeff; n++)
    {
        points[n].y = double(n);
        points[n].x = RS[n];
    }
    if (MapFuncs == 0) MapFuncs = new NaturalSpline(points, NCoeff);
    else MapFuncs->setPoints(points, NCoeff);
    nI = NI;
    nJ = NJ;
    nv = Nv;
    nChan = NChan;
    nCoeff = NCoeff;
    States = St;
    SC = C;
    Iso = I;
    R = RS;
    Data = data;
    MC = MCo;
    E = En;
    Changed();
    setImported();
}

bool CoupledSineWaveFunc::writeData(QString FN)
{
    int I, J, v, c, n;
    QFile Datei(FN);
    if (!write(&Datei)) return false;
    QTextStream S(&Datei);
    S << "Name: " << getName() << "\n";
    S << "Source: " << Source << "\n";
    for (n=0; n < nChan; n++) S << States[n]->getName() << " " << QString::number(SC[n]) << (n < nChan - 1 ? " " : "\n");
    for (n=0; n < nI; n++) S << QString::number(Iso[n]) << (n < nI - 1 ? " " : "\n");
    for (n=0; n < nCoeff; n++) S << QString::number(R[n], 'f', 7) << (n < nCoeff - 1 ? " " : "\n");
    for (I=n=0; I < nI; I++) for (J=0; J < nJ; J++) for (v=0; v < nv; v++) if (Data[I][J][v] != 0) n++;
    S << QString::number(nJ) << " " << QString::number(nv) << " " << QString::number(n) << "\n";
    for (I=0; I < nI; I++) for (J=0; J < nJ; J++) for (v=0; v < nv; v++) if (Data[I][J][v] != 0)
    {
        S << QString::number(I) << " " << QString::number(J) << " " << QString::number(v) << " " 
            << QString::number(E[I][J][v], 'f', 4) << "\n";
        for (c=0; c < nChan; c++) S << QString::number(MC[I][J][v][c], 'f', 5) << (c+1 < nChan ? " " : "\n");
        for (n=0; n < nCoeff; n++) for (c=0; c < nChan; c++) 
            S << QString::number(Data[I][J][v][c][n], 'E', 7) << (c+1 < nChan ? " " : "\n");
    }
    Saved();
    return true;
}
