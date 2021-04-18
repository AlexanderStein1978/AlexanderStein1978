//
// C++ Implementation: potential
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#include "potential.h"
#include "tools.h"
#include "MainWindow.h"
#include "termtable.h"
#include "linetable.h"
#include "molecule.h"
#include "constants.h"
#include "fit.h"
#include "inputdialog.h"
#include "fitdata.h"
#include "naturalspline.h"
#include "utils.h"
#include "duntable.h"
#include "fitobject.h"
#include "montecarlosim.h"
#include "SplinePoint.h"
#include "isotab.h"
#include "TangToenniesPot.h"
#include "AnaPot.h"
#include "PotFitLogWindow.h"
#include "tableline.h"
#include "termenergy.h"
#include "CoupledSineWaveFunc.h"
#include "splinepot.h"
#include "potentialdata.h"
#include "FitAnaPotDialog.h"
#include "CouplingFuncs.h"
#include "MLRPot.h"
#include "MTTPot.h"
#include "CoefficientDialog.h"

#include <math.h>
#include <limits>
#include <stdlib.h>

#include <QObject>
#include <QFileDialog>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QFile>
#include <QMessageBox>
#include <QDialog>
#include <QProcess>
#include <QPushButton>
#include <QLineEdit>
#include <QListWidget>
#include <QLabel>
#include <QGridLayout>
#include <QPainter>
#include <QRadioButton>
#include <QCheckBox>
#include <QMutex>
#include <QTextEdit>
#include <QPlainTextEdit>
#include <QProgressBar>


using std::numeric_limits;


Potential::Potential(MainWindow *Main, Molecule *Mol, int ThreadNum) : TableWindow(MDIChild::PotData, Main, Mol), m_saving(false),
    mWasMoving(false)
{
    setFilter("Potentials (*.pot)");
    setFileExt(".pot");
    Fit = new PotFit();
    Worker = new PotWorker(Fit);
    Fit->setWorker(Worker);
    NFC = 1;
    FixPix = new QPixmap(10, 10);
    QPainter P(FixPix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "F");
    State = 0;
    setMolecule(Mol);
    maxSplinePoints = 30;
    LRCTabOffs = adCorrTabOffs = ExcIntATabPos = -1;
    BadList = 0;
    fitN_MI = 20;
    fitN_N = 0;
    fitN_V = 4.0;
    fitN_I = 0;
    fitN_C = 0;
    fitN_A = 0.001;
    fitN_maxr = 800.0;
    fitData = 0;
    WaveFuncs = 0;
    CoupledPot = 0;
    CoupledComp = 0;
    NCoupled = 0;
    getTexTableSplineFitRunning = false;
    threadNum = ThreadNum;
    LogWindow = 0;
    m_couplingFuncs = 0;
    writePotFitTraceTab = false;
    qRegisterMetaType<FitResult>();
    Tab->setColumnCount(4);
    Tab->setHorizontalHeaderLabels(QStringList() << "coefficient" << "value" << "error" 
                                        << "sig. digits");
    Tab->setColumnWidth(0, 130);
    Tab->setColumnWidth(1, 200);
    disconnect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
    connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(updatePot()));
    connect(Fit, SIGNAL(finished()), this, SLOT(fitFinished()));
    connect(Worker, SIGNAL(badListChanged(double,double,double)), this, SLOT(badListChanged(double,double,double)));
    connect(Worker, SIGNAL(potentialImproved(PotentialData*)), this, SLOT(potentialImproved(PotentialData*)));
    /*DebugFile = new QFile(getFileName() + "DebugOut.dat");
    DebugFile->open(QIODevice::WriteOnly);
    DebugStream = new QTextStream(DebugFile);*/
}

Potential::Potential(const Potential &C) : TableWindow(PotData, C.MW, C.molecule), m_saving(false), mWasMoving(C.mWasMoving)
{
    //printf("Potential::Potential(const Potential &C)\n");
    int n;
    setFilter("Potentials (*.pot)");
    setFileExt(".pot");
    FixPix = new QPixmap(10, 10);
    QPainter P(FixPix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "F");
    Fit = new PotFit();
    Worker = C.Worker->getCopy();
    Worker->setFit(Fit);
    Fit->setWorker(Worker);
    State = C.State;
    NFC = C.NFC;
    BadList = 0;
    setName(C.getName());
    setFileName(C.getFileName());
    setMolecule(C.molecule);
    //printf("Nach setMolecule\n");
    fitN_MI = 20;
    fitN_N = 0;
    fitN_V = 4.0;
    fitN_I = 0;
    fitN_C = 0;
    fitN_A = 0.001;
    fitN_maxr = 800.0;
    getTexTableSplineFitRunning = false;
    threadNum = 0;
    LogWindow = 0;
    writePotFitTraceTab = false;
    LRCTabOffs = adCorrTabOffs = ExcIntATabPos = -1;
    /*if (C.fitData != 0)
    {
        TableLine *fD;
        fitData = new FitData(State, C.MW, C.molecule);
        fitData->setvMax(C.fitData->getvMax());
        fitData->setJMax(C.fitData->getJMax());
        C.fitData->getData(fD, n);
        fitData->setData(fD, n);
    }
    else fitData = 0;*/
    fitData = C.fitData;
    if ((WaveFuncs = C.WaveFuncs) != 0) WaveFuncs->Assign();
    if ((NCoupled = C.NCoupled) > 0)
    {
        CoupledPot = new Potential*[NCoupled];
        CoupledComp = new int[NCoupled];
        for (n=0; n < NCoupled; n++)
        {
            CoupledPot[n] = C.CoupledPot[n];
            CoupledComp[n] = C.CoupledComp[n];
        }
    }
    else
    {
        CoupledPot = 0;
        CoupledComp = 0;
    }
    m_couplingFuncs = C.m_couplingFuncs;
    Tab->setColumnCount(4);
    Tab->setHorizontalHeaderLabels(QStringList() << "coefficient" << "value" << "error" 
                                        << "sig. digits");
    Tab->setColumnWidth(0, 130);
    Tab->setColumnWidth(1, 200);
    //printf("Vor UpdateTab\n");
    UpdateTab();
    disconnect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
    connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(updatePot()));
    connect(Fit, SIGNAL(finished()), this, SLOT(fitFinished()));
    connect(Worker, SIGNAL(badListChanged(double,double,double)), this, SLOT(badListChanged(double,double,double)));
    connect(Worker, SIGNAL(potentialImproved(PotentialData*)), this, SLOT(potentialImproved(PotentialData*)));
    connect(Fit, SIGNAL(terminated()), this, SLOT(fitTerminated()));
    //printf("Ende Potential\n");
    
    /*DebugFile = new QFile("Pot" + getName() + "Debug.dat");
    DebugFile->open(QIODevice::WriteOnly);
    DebugStream = new QTextStream(DebugFile);*/
}

Potential::~Potential()
{
    delete FixPix;
    if (NCoupled > 0)
    {
        delete[] CoupledPot;
        if (CoupledComp != 0) delete[] CoupledComp;
        if (WaveFuncs != 0 ? WaveFuncs->Delete() : false) delete WaveFuncs;
    }
    delete Fit;
    delete Worker;
    //printf("Ende ~Potential\n");
    
    /*DebugFile->close();
    delete DebugStream;
    delete DebugFile;*/
}

void Potential::setWorker(PotWorker* NWorker)
{
    if (Fit->isRunning()) return;
    double *MU, *IsoF;
    IsoTab *IsoT;
    disconnect(Worker, SIGNAL(badListChanged(double,double,double)), this, SLOT(badListChanged(double,double,double)));
    disconnect(Worker, SIGNAL(potentialImproved(PotentialData*)), this, SLOT(potentialImproved(PotentialData*)));
    Worker->takeMolData(IsoT, MU, IsoF);
    delete Worker;
    Worker = NWorker;
    if (IsoT != 0) Worker->setMolData(IsoT, MU, IsoF);
    connect(Worker, SIGNAL(badListChanged(double,double,double)), this, SLOT(badListChanged(double,double,double)));
    connect(Worker, SIGNAL(potentialImproved(PotentialData*)), this, SLOT(potentialImproved(PotentialData*)));
    Fit->setWorker(NWorker);
    Worker->setFit(Fit);
}

void Potential::cutRows(int& numRows, int& numColumns, QString**& Data)
{
    setBlockChangeSignal(true);
    TableWindow::cutRows(numRows, numColumns, Data);
    setBlockChangeSignal(false);
    updatePot();
}

void Potential::DeleteRows()
{
    setBlockChangeSignal(true);
    TableWindow::DeleteRows();
    setBlockChangeSignal(false);
    updatePot();
}

void Potential::insertRows(int numRows, int numColumns, QString** Data)
{
    setBlockChangeSignal(true);
    TableWindow::insertRows(numRows, numColumns, Data);
    setBlockChangeSignal(false);
    updatePot();
}

Potential* Potential::scalePotential(double newRe, double newDe)
{
    //printf("Potential::scalePotential()\n");
    int n;
    double ARe, ADe;
    getMinimum(ARe, ADe);
    double Uinf = Worker->getUinf(), xSF = newRe / ARe, ySF = newDe / (Uinf - ADe);
    double yssSF = ySF / (xSF * xSF);
    Potential *RPot = MW->CreatePotential();
    Uinf *= ySF;
    if (Worker->isSplinePot())
    {
        SplinePoint *points, *SPoints;
        double *LRC, iA, iO, iExp, *SLRC = 0;
        int numSplinePoints, NLRC, *pLRC, *SpLRC = 0;
        Worker->getSplinePotForReading(numSplinePoints, points, NLRC, pLRC, LRC, iA, iO, iExp);
        SPoints = new SplinePoint[numSplinePoints];
        for (n=0; n < numSplinePoints; n++) 
        {
            SPoints[n].x = xSF * points[n].x;
            SPoints[n].y = ySF * points[n].y;
            SPoints[n].yss = yssSF * points[n].yss;
        }
        if (NLRC > 0)
        {
            SpLRC = new int[NLRC];
            SLRC = new double[NLRC];
            for (n=0; n < NLRC; n++) SLRC[n] = LRC[n] * pow(xSF, SpLRC[n] = pLRC[n]) * ySF;
        }
        Worker->unlock();
        iO *= ySF;
        iA *= ySF * pow(xSF, iExp) * iA;
        RPot->Worker->setSplinePot(numSplinePoints, SPoints, NLRC, SpLRC, SLRC, iA, iO, iExp);
    }
    if (Worker->isAnaPotAv() || Worker->isLPot() || Worker->isMTTPot()) 
    {
        RPot->setWorker(Worker->scalePotential(newRe, newRe));
        RPot->Worker->setFit(RPot->Fit);
    }
    RPot->Worker->setUinf(Uinf);
    RPot->setName("Scaled" + getName());
    RPot->setSource("Scaled " + getSource());
    RPot->UpdateTab();
    RPot->Changed();
    return RPot;
}

FitData* Potential::getFitData()
{
    fitData = State->getFitData();
    if (fitData == 0)
    {
        QString F = getFileName();
        QFile Fi(F = F.left(F.lastIndexOf('.')) + ".fdat");
        if(Fi.exists()) fitData = MW->getFitData(F.left(F.lastIndexOf('.')) + ".fdat", molecule);
    }
    if (fitData != 0) connect(fitData, SIGNAL(propertiesChanged()), this, SLOT(Changed()));
    return fitData;
}

void Potential::showFitData()
{
    if (fitData == 0)
    {
        QString F = getFileName();
        QFile Fi(F = F.left(F.lastIndexOf('.')) + ".fdat");
        if(Fi.exists()) fitData = MW->getFitData(F, molecule);
        if (fitData == 0)
        {
            fitData = State->getFitData();
            if (fitData == 0)
            {
                fitData = MW->CreateFitData();
                fitData->setName(getName() + "FitData");
                fitData->setFileName(F);
                State->addFitData(fitData);
            }
        }
    }
    if (!fitData->isDataAvailable()) 
    {
        fitData->updateData();
        if (!fitData->isDataAvailable())
        {
            fitData = 0;
            return;
        }
        Changed();
    }
    connect(fitData, SIGNAL(propertiesChanged()), this, SLOT(Changed()));
    if (!fitData->isVisible()) fitData->show();
    fitData->setWindowTitle("Fit data of " + getTypeString() + " " + getName());
}

void Potential::updateFitData()
{
    if (fitData == 0)
    {
        fitData = new FitData(State, MW, molecule);
        QString F = getFileName();
        QFile Fi(F = F.left(F.lastIndexOf('.')) + ".fdat");
        if (Fi.exists()) fitData->readData(F);
        connect(fitData, SIGNAL(propertiesChanged()), this, SLOT(Changed()));
        fitData->setWindowTitle("Fit data of " + getTypeString() + " " + getName());
        MW->showMDIChild(fitData);
    }
    fitData->updateData();
}

bool Potential::setElState(ElState *S)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
                                "The electronic state cannot be changed because there is currently a potential fit running!");
        return false;
    }
    State = S;
    Worker->setOmega(S->getOmega(), S->getS());
    return true;
}

void Potential::calcTermEnergies(TermTable *&TT, int NumWFPoints, int Mv, int MJ, bool show)
{
    //printf("CalcTermEnergies\n");
    int Nv, NJ, Nc, NI;
    double ****E;
    int *IsoZ = 0, *CompZ = 0;
    if (molecule == 0)
    {
        printf("Potential::calcTermEnergies: Error: no molecule available!\n");
        return;
    } 
    if (MW == 0)    
    {
        printf("Potential::calcTermEnergies: Error: no MainWindow available!\n");
        return;
    }
    if (WaveFuncs != 0)
    {
        //int C, L = State->getLambda(), n, m;
        //ElState *St;
        int *Z, n;
        E = WaveFuncs->getEnergy();
        WaveFuncs->getIso(NI, Z);
        IsoZ = new int[NI];
        for (n=0; n < NI; n++) IsoZ[n] = Z[n];
        CompZ = new int[Nc = 1];
        CompZ[0] = 0;
        Nv = WaveFuncs->getNv();
        NJ = WaveFuncs->getNJ();
        /*if (L==0)
        {
            for (n=m=0; n < Nc; n++)
            {
                WaveFuncs->getStateComp(n, St, C);
                if (St->getLambda() == 0 || C%2 == 0) 
                {
                    CompZ[m] = m;
                    m++;
                }
                else
                {
                    Destroy(E[n], NI, Nv);
                    E[n] = 0;
                }
            }
            if (m<n)
            {
                for (n=0; E[n] != 0; n++) ;
                for (m=n+1; m < Nc; m++) if (E[m] != 0) E[n++] = E[m];
                Nc = n;
            }
        }
        else for (n=m=0; n < Nc; n++)
        {
            WaveFuncs->getStateComp(n, St, C);
            if (St->getLambda() == 0 ? m%2 == 1 : (m*C)%2 == 1) m++; 
            CompZ[n] = m++;
        }*/
    }
    else
    {
        Nc = ((State != 0 ? State->getLambda() : 1) > 0 ? 2 : 1);
        if (Worker->getNSpinRGamma() > 0) Nc *= (int(2.0 * State->getS()) + 1);
        E = new double***[Nc];
        CompZ = new int[Nc];
        TermTable *cTT = (State != 0 ? State->getTermTable(false) : 0);
        int n;
        if (cTT != 0 ? cTT->getNumComp() == Nc : false)
        {
            int *cCompZ = cTT->getCompZ();
            for (n=0; n < Nc; n++) CompZ[n] = cCompZ[n];
        }
        else for (n=0; n < Nc; n++) CompZ[n] = n;
        if (Mv < 0 || MJ < 0)
        {
            InputDialog Dialog(70, 350); 
            Dialog.setLabel("Maximum v", "Maximum J");
            if (Dialog.exec() == QDialog::Rejected) return;
            Dialog.getData(Mv, MJ);
        }
        getEV(E, Nc, NI = Worker->getNIso(), Nv = Mv + 1, NJ = MJ + 1, NumWFPoints);
    }
    if (E == 0)
    {
        delete[] E;
        return;
    }
    if (TT == 0)
    {
        TT = MW->CreateTermTable();
        if (State != 0) State->addTermTable(TT);
        TT->setName("Term" + getName());
        TT->setFileName("Term" + getFileName());
        //printf("Vor setSource\n");
        TT->setSource("Calculated from " + getSource());
    }
    //printf("Nach getEV\n");
    TT->setData(E, Nc, NI, Nv - 1, NJ - 1, CompZ);
    if (IsoZ != 0) TT->setIsoZ(IsoZ);
    if (show) TT->show();
    //printf("Ende CalcTermEnergies, Nv=%d\n", Nv);
}

void Potential::getFCF(Potential *uPot, int Iso, int Js, int Jss, int &Nvs, int &Nvss, double ***&FCF, int NumWFPoints)
{
    //printf("Potential::getFCF, Iso=%d, Js=%d, Jss=%d\n", Iso, Js, Jss);
    int *Pu, *Px, *rix, *rax, *riu, *rau, UWI = -1, n;
    double *Eu, *Ex, **WFu, **WFx, *Fu, *Fx, *Nu, *Nx, Min = getMinimumE(), uMin = uPot->getMinimumE();
    CoupledSineWaveFunc *UWF = uPot->getCoupledSineWaveFunc();
    getWaveFunc(Iso, Jss, Nvss, Ex, WFx, rix, rax, Px, Fx, Nx, NumWFPoints);
    if (UWF != 0 ? UWF->getNJ() > Js : false)
    {
        int N, *I;
        UWF->getIso(N, I);
        for (n=0; n<N; n++) if (I[n] == Iso) UWI = n;
    }
    if (UWI >= 0)
    {
        int N, *v, c;
        ElState *Stu = uPot->getElState(), *Stc = 0;
        UWF->getv(UWI, Js, N, v);
        Pu = new int[Nvs = v[N-1] + 1];
        riu = new int[Nvs];
        rau = new int[Nvs];
        Eu = new double[Nvs];
        WFu = new double*[Nvs];
        Fu = new double[Nvs];
        Nu = new double[Nvs];
        for (n=0, N = UWF->getNChannels(); Stc != Stu && n<N; n++) UWF->getStateComp(n, Stc, c);
        c=n-1;
        for (n=0; n < Nvs; n++)
        {
            riu[n] = Pu[n] = 0;
            rau[n] = NumWFPoints - 1;
            Fu[n] = Nu[n] = 1.0;
            WFu[n] = new double[NumWFPoints];
        }
        if (!UWF->getWaveFuncs(UWI, Js, c, rmin, rmax, NumWFPoints, WFu, Eu)) for (n=0; n < Nvs; n++)
        {
            for (N=0; N < NumWFPoints; N++) WFu[n][N] = 0.0;
            Eu[n] = 0.0;
        }
    }
    else uPot->getWaveFunc(Iso, Js, Nvs, Eu, WFu, riu, rau, Pu, Fu, Nu, NumWFPoints);
    //printf("getFCF: Nvss=%d\n", Nvss);
    //printf("Vor Create FCF\n");
    getFCF(Nvss, Nvs, Ex, Eu, WFx, WFu, rix, riu, rax, rau, Px, Pu, Fx, Fu, Nx, Nu, FCF, Min, uMin);
    //printf("Name=%s, uName=%s\n", getName().ascii(), uPot->getName().ascii());
}

void Potential::getFCF(int &Nvss, int &Nvs, double* Ex, double* Eu, double** WFx, double** WFu, 
                       int* rix, int* riu, int* rax, int* rau, int* Px, int* Pu, double* Fx, 
                       double* Fu, double* Nx, double* Nu, double ***&FCF, double Min, double uMin)
{
    //printf("Ex[1]=%f, Eu[1]=%f\n", Ex[1], Eu[1]);
    int Mvs = Nvs, Mvss = Nvss, i, vs, vss, P1, P2, ra;
    double GS, S, F1, F2;
    //printf("Min=%f, uMin=%f\n", Min, uMin);
    FCF = Create(Nvs, Nvss, 2);
    //printf("Nach Create FCF\n");
    for (vss = 0; vss < Nvss; vss++) for (vs = 0; vs < Nvs; vs++)
    {
        if (Pu[vs] < Px[vss])
        {
            P1 = Pu[vs];
            P2 = Px[vss];
            F1 = Fu[vs];
            F2 = F1 * Fx[vss];
        }
        else
        {
            P1 = Px[vss];
            P2 = Pu[vs];
            F1 = Fx[vss];
            F2 = F1 * Fu[vs];
        }
        if ((ra = (rax[vss] < rau[vs] ? rax[vss] : rau[vs])) < P1) P1 = ra - 1;
        if (ra < P2) P2 = ra - 1;
        for (i = (rix[vss] > riu[vs] ? rix[vss] : riu[vs]), GS = 0.0; i <= P1; i++) 
            GS += WFu[vs][i] * WFx[vss][i];
        for (S=0.0; i <= P2; i++) S += WFu[vs][i] * WFx[vss][i];
        GS += S * F1;
        for (S=0.0; i < ra; i++) S += WFu[vs][i] * WFx[vss][i];
        GS += S * F2;
        //if (vs == 250 && vss == 60) 
            //printf("ri=%d, ra=%d, P1=%d, P2=%d, GS=%f, Nu=%f, Nx=%f\n", ri, ra, P1);
        FCF[vs][vss][0] = GS * GS * Nu[vs] * Nx[vss];
        FCF[vs][vss][1] = Eu[vs] - Ex[vss];
        if (Min < 0.0) FCF[vs][vss][1] += Min;
        if (uMin < 0.0) FCF[vs][vss][1] -= uMin;
    }
    //printf("Vor delete E\n");
    delete[] Eu;
    delete[] Ex;
    delete[] rix;
    delete[] rax;
    delete[] riu;
    delete[] rau;
    //printf("Vor delete P\n");
    delete[] Pu;
    delete[] Px;
    //printf("Vor delete F\n");
    delete[] Fu;
    delete[] Fx;
    //printf("Vor delete N\n");
    delete[] Nu;
    delete[] Nx;
    //printf("Vor Destroy WF\n");
    Destroy(WFu, Mvs);
    Destroy(WFx, Mvss);
    //printf("Ende GetFCF, Nvs=%d\n", Nvs);
}

double *Potential::calcSigFuncs(int /*numPoints*/, double /*rMin*/, double /*rMax*/)
{
    return 0;
}

void Potential::getWaveFunc(int I, int J, int &Nv, double *&EL, double **&WF, int *&ri, int *&ra, 
                            int *&CP, double *&F, double *&SQ, int NumWFPoints)
{
    //printf("Potential::getWaveFunc, Nv=%d\n", Nv);
    double h = (rmax - rmin) / (NumWFPoints - 1), *RPot = new double[NumWFPoints], *Y = new double[NumWFPoints];
    double *U = getPoints(rmin, rmax, NumWFPoints), *IsoF = Worker->getIsoF_ForReading(), MF = IsoF[I] * h*h;
    int n, nv = Nv;
    WF = Create(Nv, NumWFPoints);
    //printf("Nach Create WF\n");
    //aPot->GetMin(R, M);
    F = new double[Nv];
    SQ = new double[Nv];
    EL = new double[Nv];
    CP = new int[Nv];
    ri = new int[Nv];
    ra = new int[Nv];
    Worker->CRotPot(J, I, U, RPot, rmin, h, NumWFPoints);
    NumerovCooley(Y, MF, Nv, EL, WF, ri, ra, CP, F, SQ, NumWFPoints, RPot);
    delete[] U;
    delete[] Y;
    delete[] RPot;
    for (n = Nv; n < nv; n++) delete[] WF[n];
    //printf("Ende von getWaveFunc, CP[0]=%d, Nv=%d\n", CP[0], Nv);
}

void Potential::exportWaveFunction(int NumWFPoints)
{
    if (WaveFuncs == 0) return;
    int n, NI, *I, NC = WaveFuncs->getNChannels(), v, J, i, Nv = WaveFuncs->getNv();
    double ***WF = Create(NC, Nv, NumWFPoints), E[Nv], r, h = (rmax - rmin) / ((NumWFPoints - 1) * a0_Angstrom);
    double UF = 0.7270458;
    IsoTab *Iso = molecule->getIso();
    QDialog *D = new QDialog(this);
    QGridLayout *L = new QGridLayout(D);
    QComboBox *IB = new QComboBox(D);
    QLineEdit *vE = new QLineEdit("0", D), *JE = new QLineEdit("0", D);
    QPushButton *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    QString FileName;
    D->setWindowTitle("MolSpektAnalysis");
    L->addWidget(new QLabel("Please select the level of which the wave function shall be exported:", D), 0, 0, 1, 2);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Iso:", D), 2, 0);
    L->addWidget(IB, 2, 1);
    L->addWidget(new QLabel("v:", D), 3, 0);
    L->addWidget(vE, 3, 1);
    L->addWidget(new QLabel("J:", D), 4, 0);
    L->addWidget(JE, 4, 1);
    L->setRowMinimumHeight(5, 20);
    L->addWidget(OK, 6, 0);
    L->addWidget(Cancel, 6, 1);
    WaveFuncs->getIso(NI, I);
    for (n=0; n < NI; n++) IB->addItem(Iso->getIsoName(I[n]));
    IB->setEditable(false);
    connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
    if (D->exec() == QDialog::Accepted)
    {
        i = IB->currentIndex();
        v = vE->text().toInt();
        J = JE->text().toInt();
        for (n=0; n < Nv; n++) E[n] = 0.0;
        if (v >= 0 && v < Nv && J >= 0 ? WaveFuncs->getWaveFuncs(i, J, rmin, rmax, NumWFPoints, WF, E) : false)
        {
            if (E[v] != 0.0)
            {
                if (!(FileName = QFileDialog::getSaveFileName(this, "Select file name", MW->getDir(), 
                                                                    "Text files (*.dat)")).isEmpty()) 
                {
                    QFile F(FileName);
                    F.open(QIODevice::WriteOnly);
                    QTextStream S(&F);
                    S << "Isotopologue " << Iso->getIsoName(I[i]) << ", v'=" << QString::number(v) << ", J'=" 
                      << QString::number(J) << "\n"; 
                    for (n=0, r = rmin / a0_Angstrom; n < NumWFPoints; n++, r+=h)
                    {
                        S << QString::number(r, 'f', 12);
                        for (i=0; i < NC; i++) S << "\t" << QString::number(WF[i][v][n] * UF, 'e', 12);
                        S << "\n";
                    }
                }
            }
            else QMessageBox::information(this, "MolSpektAnalysis", "Error: the level with v'=" + QString::number(v) +
                                            + "could not be found!");
        }
        else QMessageBox::information(this, "MolSpektAnalysis", "Error: the level with v'=" + QString::number(v)
                                        + " and J'=" + QString::number(J) + " does not exist!");
    }
    Destroy(WF, NC, Nv);
    delete Iso;
    delete D;
}

void Potential::getExchangeInt(double& A, double& alpha, double& gamma)
{
    AnaPot *aPot = Worker->getAnaPotForReading();
    if (aPot != 0) 
    {
        aPot->getExchangeCoeff(A, alpha, gamma);
        Worker->unlock();
    }
    else A = alpha = gamma = 0.0;
}

void Potential::setExchangeInt(double A, double alpha, double gamma)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
            "The spin rotation coefficients cannot be changed because there is currently a fit of the potential running!");
        return;
    }
    AnaPot *aPot = Worker->getAnaPotForWriting();
    if (aPot != 0) aPot->setExchangeCoeff(A, alpha, gamma);
}

void Potential::getAdCorrForWriting(int& NAdCorr, double*& adCorr, int& TAdCorr, int& PAdCorr, double& RIso1, double& RIso2, 
                                    double& Rm, double& b)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
                        "The adiabatic correction cannot be modified because there currentliy is a fit of the potential running!");
        NAdCorr = 0;
    }
    else Worker->getAdCorrForWriting(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b);
}

bool Potential::setAdCorr(int NAdCorr, double* adCorr, int TAdCorr, int PAdCorr, double RIso1, double RIso2, double Rm, double b)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
                        "The adiabatic correction cannot be modified because there currentliy is a fit of the potential running!");
        return false;
    }
    return Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b);
}

void Potential::getLRCoeffForReading(int& N, int*& pLRC, double*& RLRC)
{
    Worker->getLRCoeffForReading(N, pLRC, RLRC);
}

void Potential::getLRCoeffForWriting(int& N, int*& pLRC, double*& RLRC, bool **LRCFree)
{
    N=0;
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
            "The long range coefficients cannot be changed because there is currently a fit of the potential running!");
        return;
    }
    Worker->getLRCoeffForWriting(N, pLRC, RLRC, LRCFree);
    
}

void Potential::getSplinePotForWriting(int &NSplinePoints, SplinePoint *&points, int &NLRC, int *&pLRC, double *&LRC, double &iA,
    double &iO, double &Exp)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
            "The potential cannot be modified because there is currently a fit of the potential running!");
        NSplinePoints = 0;
        points = nullptr;
        NLRC = 0;
        pLRC = nullptr;
        LRC = nullptr;
        iA = iO = Exp = 0;
        return;
    }
    Worker->getSplinePotForWriting(NSplinePoints, points, NLRC, pLRC, LRC, iA, iO, Exp);
}

void Potential::getFastWaveFunc(int I, int J, int C, double E, double M, double h, int NP, double *WF, 
                                double &Ri, double &Ra, double &ISq, double &F, int &i, int &c, int &a)
{
    //printf("Potential::getFastWaveFunc\n");
    double *U = new double[NP], *P = getPoints(Ri, Ra, NP), PM = M * 0.0833333333333, *A = new double[NP];
    if (Worker->isAdCorrA()) 
    {
        IsoTab *IsoT = Worker->getIsoT();
        Worker->CalcAdCorrPot(P, A, Ri, h, NP, IsoT->mIso1[I], IsoT->mIso2[I]);
        delete[] P;
    }
    else A=P;
    Worker->CRotPot(J, I, A, U, Ri, h, NP, C);
    delete[] A;
    Worker->fastWaveFunc(U, E, PM, h, NP, WF, Ri, Ra, ISq, F, i, c, a);
    delete[] U;
}

bool Potential::getFastWaveFuncs(int C, int I, int J, int Nv, double M, double h, int NP, double** WF, 
                                 double ri, double ra, double* WFsq, double* F, int* i, int* cp, int* a)
{
    TermTable *TT = State->getTermTable();
    if (C >= TT->getNumComp() || I >= TT->getNumIso() || J > TT->getMaxJ() 
        || Nv > TT->getMaxv() + 1) return false;
    double ****EData = TT->getData(), *U = getPoints(ri, ra, NP, C), *A, *R, PM = M * 0.0833333333333, Ri, Ra;
    int fc, c, l, v;
    float S;
    IsoTab *IsoT = Worker->getIsoT();
    if ((S = State->getS()) > 0.0 && (l = State->getLambda()) > 0)
    {
        fc = C/2;
        c = C - 2 * fc;
    }
    else if (S > 0.0)
    {
        fc = C;
        c=0;
    }
    else
    {
        c = C;
        fc = 0;
    }
    if (Worker->isAdCorrA())
    {
        Worker->CalcAdCorrPot(U, A = new double[NP], ri, h, NP, IsoT->mIso1[I], IsoT->mIso2[I]);
        delete[] U;
    }
    else A=U;
    Worker->CRotPot(J, I, A, R = new double[NP], ri, h, NP, c, fc);
    delete[] A;
    for (v=0; v < Nv; v++) 
        Worker->fastWaveFunc(R, EData[C][I][v][J], PM, h, NP, WF[v], Ri = ri, Ra, WFsq[v], F[v], i[v], cp[v], a[v]);
    delete[] R;
    return true;
}

void Potential::getFastFCF(int I, Potential *uPot, int Js, int Comp, double uE, int N, double *lE, int *lJss, 
                           double *FCF, int NumWFPoints, double **uWF, double *TStr, int NWFC)
{
    //printf("Potential::getFastFCF, N=%d, uE=%f\n", N, uE);
    if (NWFC < 1) NWFC = 1;
    int NP = NumWFPoints, i, j, n, ui, ua, li, la, s, uc, c, UWI = -1, vs, wf;
    double Ri = rmin, Ra = rmax, uISq, ISq, Sq, WF, uF, F, fcf[NWFC], fcf2[NWFC], tFCF[NWFC];
    double h = (rmax - rmin) / double(NP - 1), uRi = Ri, uRa = Ra;
    double *U = new double[NP], eins = 1.0, *IsoF = Worker->getIsoF_ForReading(), M = IsoF[I] * h*h, Y0, Y1, Y2 = 0.0;
    double PM = M * 0.0833333333333, ED, *P = getPoints(Ri, Ra, NP);
    //printf("Vor lE\n");
    if ((ED = getMinimumE()) < 0.0) for (i=0; i<N; i++) lE[i] += ED;
    
    /*int Nv = 60, *ri, *ra, *CP;
    double **tWF, *EL, *tF, *SQ;
    printf("h=%g\n", h);
    getWaveFunc(I, lJss[0], Nv, EL, tWF, ri, ra, CP, tF, SQ);
    delete[] ri;
    delete[] ra;
    delete[] CP;
    Destroy(tWF, Nv);
    delete[] tF;
    delete[] SQ;
    for (n=0, i=-1; n < Nv; n++)
    {
        printf("v=%d, Ec=%g, ", n, EL[n]);
        for (i++; (i<N ? lJss[i] != lJss[0] : false); i++) ;
        printf("Ei=%g, D=%g\n", lE[i], EL[n]-lE[i]);
    }
    delete[] EL;
    uPot->getWaveFunc(I, Js, Nv, EL, tWF, ri, ra, CP, tF, SQ);
    delete[] ri;
    delete[] ra;
    delete[] CP;
    Destroy(tWF, Nv);
    delete[] tF;
    delete[] SQ;
    for (n=0; n < Nv; n++) printf("Eu=%g, E[%d]=%g, D=%g\n", uE, n, EL[n], uE-EL[n]);
    delete[] EL;*/
    
    //printf("Vor getFastWF\n");
    if (uWF == 0)
    {
        CoupledSineWaveFunc *UWF = uPot->getCoupledSineWaveFunc();
        if (UWF != 0 ? UWF->getNJ() > Js : false)
        {
            int N, *Iso;
            UWF->getIso(N, Iso);
            for (n=0; n<N; n++) if (Iso[n] == I) UWI = n;
        }
        if (UWI >= 0 ? (vs = UWF->getv(UWI, Js, uE)) >= 0 : false)
        {
            int N, c;
            double EB;
            ElState *Stu = uPot->getElState(), *Stc = 0;
            uWF = Create(1, NP);
            for (n=0, N = UWF->getNChannels(); Stc != Stu && n<N; n++) UWF->getStateComp(n, Stc, c);
            c=n-1;
            UWF->getWaveFunc(UWI, Js, vs, c, rmin, rmax, NumWFPoints, uWF[0], EB, N);
        }
    }
    if (uWF != 0)
    {
        uRi = rmin;
        uRa = rmax;
        uISq = 1.0;
        uF = 1.0;
        ui = uc = 0;
        ua = NumWFPoints - 1;
    }
    else 
    {
        uWF = Create(1, NP);
        uPot->getFastWaveFunc(I, Js, Comp, uE, M, h, NP, uWF[0], uRi, uRa, uISq, uF, ui, uc, ua);
    }
    //printf("Nach getFastWF\n");
    for (i=0; i<N; i++) FCF[i] = -1.0;
    for (i=0; i<N; i++) if (FCF[i] == -1.0)
    {
        Worker->CRotPot(lJss[i], I, P, U, Ri, h, NP);
        for (j=i; j<N; j++) if (lJss[j] == lJss[i])
        {
            //printf("j=%d\n", j);
            getMinMaxInt(U, lE[j], M, NP, li, la);
            //printf("Nach gerMinMaxInt\n");
            if (la < NP - 1 ? U[la] < lE[j] : true)
            {
                for (la = NP - 1; (la > li? U[la - 1] > U[la] || U[la] < lE[j] : false); la--) ;
                if (la == li) 
                {
                    printf("E[%d]=%f is too high!\n", j, lE[j]);
                    FCF[j] = 0.0;
                    continue;
                }
            }
            Y1 = 1e-5;
            Y0 = Y1 * exp(-sqrt(M*(U[la] - lE[j])));
            WF = Y0 / (eins - PM * (U[c = la] - lE[j]));
            Sq = WF * WF;
            if (la < ua) for (wf = 0; wf < NWFC; wf++) 
            {
                fcf[wf] = WF * uWF[wf][c];
                if (la <= uc) fcf[wf] /= uF;
            }
            else for (wf = 0; wf < NWFC; wf++) fcf[wf] = 0.0;
            WF = Y1 / (eins - PM * (ED = U[--c] - lE[j]));
            Sq += WF * WF;
            if (c < ua) for (wf = 0; wf < NWFC; wf++)
            {
                if (c > uc) fcf[wf] += WF * uWF[wf][c];
                else fcf[wf] += WF * uWF[wf][c] / uF;
            }
            for (c--; c > ua; c--)
            {
                Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                WF = Y2 / (eins - PM * (ED = U[c] - lE[j]));
                if (Y2 <= Y1) break;
                Sq += WF * WF;
                Y0 = Y1;
                Y1 = Y2;
            }
            if (Y2 > Y1)
            {
                for (c--; c > uc; c--)
                {
                    Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                    WF = Y2 / (eins - PM * (ED = U[c] - lE[j]));
                    if (Y2 <= Y1) break;
                    Sq += WF * WF;
                    for (wf = 0; wf < NWFC; wf++) fcf[wf] += WF * uWF[wf][c];
                    Y0 = Y1;
                    Y1 = Y2;
                }
                for (wf = 0; wf < NWFC; wf++) fcf[wf] *= uF;
                if (Y2 > Y1)
                {
                    for (c--; c > ui; c--)
                    {
                        Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                        WF = Y2 / (eins - PM * (ED = U[c] - lE[j]));
                        if (Y2 <= Y1) break;
                        Sq += WF * WF;
                        for (wf = 0; wf < NWFC; wf++) fcf[wf] += WF * uWF[wf][c];
                        Y0 = Y1;
                        Y1 = Y2;
                    }
                    if (Y2 > Y1)
                    {
                        for (c--; c > li; c--)
                        {
                            Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                            WF = Y2 / (eins - PM * (ED = U[c] - lE[j]));
                            if (Y2 <= Y1) break;
                            Sq += WF * WF;
                            Y0 = Y1;
                            Y1 = Y2;
                        }
                    }
                }
            }
            F = eins / WF;
            //printf("ui=%d, li=%d, ua=%d, la=%d, NP=%d\n", ui, li, ua, la, NP);
            Y0 = (Y1 = 1e-5) * exp(-sqrt(M * (U[li] - lE[j])));
            WF = Y0 / (eins - PM * (U[li] - lE[j]));
            ISq = WF * WF;
            if (li >= ui) for (wf = 0; wf < NWFC; wf++)
            {
                tFCF[wf] = WF * uWF[wf][li];
                if (li > uc) tFCF[wf] *= uF;
            }
            else for (wf = 0; wf < NWFC; wf++) tFCF[wf] = 0.0;
            WF = Y1 / (eins - PM * (ED = U[n = li + 1] - lE[j]));
            if (n >= ui) for (wf = 0; wf < NWFC; wf++)
            {
                if (n > uc) tFCF[wf] += WF * uWF[wf][n] * uF;
                else tFCF[wf] += WF * uWF[wf][n];
            }
            ISq += WF * WF;
            for (n++, s = (ui <= c ? ui : c+1); n<s; n++)
            {
                Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                WF = Y2 / (eins - PM * (ED = U[n] - lE[j]));
                ISq += WF * WF;
                Y0 = Y1;
                Y1 = Y2;
            }
            for (s = (uc < c ? uc : c); n<=s; n++)
            {
                Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                WF = Y2 / (eins - PM * (ED = U[n] - lE[j]));
                ISq += WF * WF;
                for (wf = 0; wf < NWFC; wf++) tFCF[wf] += WF * uWF[wf][n];
                Y0 = Y1;
                Y1 = Y2;
            }
            for (wf = 0; wf < NWFC; wf++) fcf2[wf] = 0.0;
            for (s = (ua <= c ? ua : c+1); n<s; n++)
            {
                Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                WF = Y2 / (eins - PM * (ED = U[n] - lE[j]));
                ISq += WF * WF;
                for (wf = 0; wf < NWFC; wf++) fcf2[wf] += WF * uWF[wf][n];
                Y0 = Y1;
                Y1 = Y2;
            }
            while (n <= c)
            {
                Y2 = M * ED * WF + 2.0 * Y1 - Y0;
                WF = Y2 / (eins - PM * (ED = U[n++] - lE[j]));
                ISq += WF * WF;
                Y0 = Y1;
                Y1 = Y2;
            }
            F *= WF;
            ISq += F*F * Sq;
            for (wf = 0, FCF[j] = 0.0; wf < NWFC; wf++)
            {
                tFCF[wf] += fcf2[wf] * uF + fcf[wf] * F;
                //printf("n=%d, la=%d\n", n, la);
                //printf("FCF=%g, ISq=%g, uISq=%g\n", FCF[j], ISq, uISq);
                tFCF[wf] *= tFCF[wf];
                tFCF[wf] /= (ISq * uISq);
                FCF[j] += tFCF[wf] * (TStr !=0 ? TStr[wf] : 1.0);
            }
                //printf("lJss[%d]=%d, lE[%d]=%f, FCF[%d]=%f\n", j, lJss[j], j, lE[j], j, FCF[j]);
        }
    }
    delete[] P;
    Destroy(uWF, NWFC);
    delete[] U;
    //printf("Ende getFastFCF\n");
}

void Potential::exportAsymptoticLevels(QString FileName, int numv, int maxJ, int NumWFPoints)
{
    if (State == 0 || molecule == 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Error: the data cannot be exported!");
        return;
    }
    double *U = getPoints(rmin, rmax, NumWFPoints);
    if (U==0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Error: the data cannot be exported!");
        return;
    }
    int J, nv, I, v, Nv, NJ, NI, NIso = Worker->getNIso(), JStart[NIso], JStep[NIso], Rm, RM, LCP;
    double h = (rmax - rmin) / (NumWFPoints - 1), Max = getAsymptote(), *RPot = new double[NumWFPoints], MF = 0.0;
    double hsq = h*h, DE, *Y = new double[NumWFPoints], *IsoF = Worker->getIsoF_ForReading();
    for (I=0; I < NIso; I++) 
    {
        if (IsoF[I] > MF) MF = IsoF[I];
        JStart[I] = State->getJStart(I, 0);
        JStep[I] = molecule->getJStep(I);
    }
    MF *= hsq;
    getMinMaxInt(U, Max, MF, NumWFPoints, Rm, RM);
    double *LWF = new double[RM - Rm], LF, LSQ;
    Numerov(MF, Max, DE, Nv, RM - Rm, U + Rm, LWF, LCP, LF, LSQ, Y);
    delete[] LWF;
    if (DE > 0.0) Nv--;
    double F[Nv], SQ[Nv], **WF = Create(Nv, NumWFPoints);
    double ***EL = new double**[NI = NIso];
    int Maxv = Nv - 1, MaxJ = maxJ, CP[Nv], ri[Nv], ra[Nv];
    for (I=0, NJ = 0, Nv = 0; I < NI; I++) 
    {
        EL[I] = new double*[MaxJ+1];
        MF = IsoF[I] * hsq;
        for (J = JStart[I]; J <= MaxJ; J += JStep[I])
        {
            EL[I][J] = new double[Maxv+1];
            Worker->CRotPot(J, I, U, RPot, rmin, h, NumWFPoints);
            NumerovCooley(Y, MF, nv = Maxv + 1, EL[I][J], WF, ri, ra, CP, F, SQ, NumWFPoints, RPot);
            printf("Iso=%d, J=%d, Maxv=%d\n", I, J, nv-1);
            if (nv <= Maxv) EL[I][J][nv] = 0.0;
            if (nv > Nv) Nv = nv;
            else if (nv == 0)
            {
                delete[] EL[I][J];
                EL[I][J] = NULL;
                break;
            }
        }
        if (J > NJ) NJ = J;
    }
    if (NJ > MaxJ + 1) NJ = MaxJ + 1;
    Destroy(WF, Maxv + 1);
    delete[] U;
    delete[] RPot;
    delete[] Y;
    QFile File(FileName);
    File.open(QIODevice::WriteOnly);
    QTextStream S(&File);
    IsoTab *IsoT = Worker->getIsoT();
    for (I=0; I<NI; I++)
    {
        S << "Isotopologue " + IsoT->getIsoName(I) + ":\n";
        for (J=JStart[I]; (J < NJ ? EL[I][J] != NULL : false); J += JStep[I])
        {
            S << "J=" + QString::number(J) + ":\n";
            for (v = Nv - numv; (v < Nv ? EL[I][J][v] != 0.0 : false); v++) 
                S << "v = " + QString::number(v) + " E = " << QString::number(EL[I][J][v] - Max, 'f', 4) + "\n";
            delete[] EL[I][J];
        }
        delete[] EL[I];
        S << "\n";
    }
    delete[] EL;
}

void Potential::getEV(double ****&Ev, int &Nc, int &NI, int &Nv, int &NJ, int NumWFPoints)
{
//     printf("Beginn getEV\n");
    bool AdCorr = Worker->isAdCorrA();
    double *U = getPoints(rmin, rmax, NumWFPoints);
    if (U==0)
    {
        Ev = 0;
        NI = Nv = NJ = Nc = 0;
        return;
    }
    double h = (rmax - rmin) / (NumWFPoints - 1), *RPot = new double[NumWFPoints], **WF = Create(Nv, NumWFPoints);
    double *IsoF = Worker->getIsoF_ForReading();
    double F[Nv], SQ[Nv], hsq = h*h, MF, *AdCorrPot = (AdCorr ? new double[NumWFPoints] : U), *Y = new double[NumWFPoints];
    int nc = (Worker->isLambdaDoublingA() && Nc > 1 ? 2 : 1);
    double ****EL = new double***[Nc = ((State != 0 ? State->getS() : 0.0) > 0.0 && Nc >= 2 * nc ? 2 * nc : nc)], R, M;
    int J, nv, I, v, Maxv = Nv - 1, MaxJ = NJ - 1, CP[Nv], ri[Nv], ra[Nv], c, FC, ef;
    for (c=0, NI = Worker->getNIso(); c < Nc; c++) EL[c] = new double**[NI];
    IsoTab *IsoT = Worker->getIsoT();
    printf("Nv=%d, Maxv=%d\n", Nv, Maxv);
    //NI=1;
    getMinimum(R, M);
    /*QFile DDatei("DebugTerm.dat");
    DDatei.open(QIODevice::WriteOnly);
    QTextStream S1(&DDatei);*/
    //printf("NIso=%d, MaxJ=%d, Maxv=%d\n", NI, MaxJ, Maxv);
        
    for (c=0; c < Nc; c++) 
    {
        FC = (Nc > nc ? c / nc + 1 : 0);
        ef = (Nc > nc ? c - nc * (FC - 1) : c);
        for (I=0, NJ = 0, Nv = 0; I < NI; I++) 
        {
            //S1 << "Iso=" << I << "\n";
            EL[c][I] = new double*[MaxJ+1];
            MF = IsoF[I] * hsq;
            if (AdCorr) Worker->CalcAdCorrPot(U, AdCorrPot, rmin, h, NumWFPoints, IsoT->mIso1[I], IsoT->mIso2[I]);
            for (J=0; J <= MaxJ; J++)
            {
                //S1 << "J=" << J << "\n";
                EL[c][I][J] = new double[Maxv+1];
                Worker->CRotPot(J, I, AdCorrPot, RPot, rmin, h, NumWFPoints, ef, FC);
                //printf("Vor Numerov\n");
                NumerovCooley(Y, MF, nv = Maxv + 1, EL[c][I][J], WF, ri, ra, CP, F, SQ, NumWFPoints, RPot);
                printf("Iso=%d, J=%d, Maxv=%d\n", I, J, nv-1);
                if (nv <= Maxv) EL[c][I][J][nv] = 0.0;
                if (nv > Nv) Nv = nv;
                else if (nv == 0)
                {
                    delete[] EL[c][I][J];
                    EL[c][I][J] = 0;
                    break;
                }
                //for (v=0; v<nv; v++) 
                    //S1 << "v=" << v << "\tT=" << QString::number(EL[I][J][v] - M, 'g', 10) << "\n";
            }
            if (J > NJ) NJ = J;
        }
    }
    //printf("Nach Schleife, Nv=%d, Maxv=%d\n", Nv, Maxv);
    //DDatei.close();
    
    Destroy(WF, Maxv + 1);
    //printf("Nach Destroy WF\n");
    delete[] U;
    delete[] RPot;
    delete[] Y;
    if (AdCorr) delete[] AdCorrPot;
    //printf("Nach delete U\n");
    Ev = Create(Nc, NI, Nv, NJ);
    for (c=0; c < Nc; c++) 
    {
        for (I=0; I<NI; I++)
        {
            for (J=0; (J < NJ ? EL[c][I][J] != NULL : false); J++)
            {
                for (v=0; (v < Nv ? EL[c][I][J][v] != 0.0 : false); v++) 
                    Ev[c][I][v][J] = (M>0 ? EL[c][I][J][v] : EL[c][I][J][v] - M);
                for (; v < Nv; v++) Ev[c][I][v][J] = 0.0;
                //printf("Vor delete EL[%d][%d]\n", I, J);
                delete[] EL[c][I][J];
            }
            for (; J < NJ; J++) for (v=0; v < Nv; v++) Ev[c][I][v][J] = 0.0;
            //printf("Vor delete EL[%d]\n", I);
            delete[] EL[c][I];
        }
        delete[] EL[c];
    }
    //printf("Vor delete EL\n");
    delete[] EL;
    //printf("Ende von getEv, Nv=%d, Min=%f, E[9][1][60]=%f\n", Nv, M, Ev[9][1][60]);
}

void Potential::setFitData(int maxv, int maxJ, int Iso, bool OnlyTerm)
{
    if (Fit->isRunning()) return;
    int *mJ, **mv, Mv, nES, *nEL, nGS, *nGP, *nGL, **PS, MaxJ, Maxv, NIso = Worker->getNIso(); 
    TermEnergy **EL;
    TableLine **GL; 
    nES = nGS = 0;
    int i, e, g, n, m, j, N, DL[MaxLineTables];
    TableLine *Lines;
    LineTable *LineT = 0;
    FitData *fDB = (State != 0 ? State->getFitData() : 0);
    if (fDB != 0) fitData = fDB;
    if (fitData == 0) 
    {
        showFitData();
        if (fitData == 0) return;
    }
    fitData->getData(Lines, N, -1, -2, -2, 0, Iso, State);
    if (N == 0) return;
    for (i=0; i < MaxLineTables; i++) DL[i] = 0;
     maxJ = MaxJ = fitData->getJMax();
    maxv = Maxv = fitData->getvMax();
    for (n=m=0; n<N; n++) if (Lines[n].LTab != LineT || (n>0 ? Lines[n-1].isTE != Lines[n].isTE : true))
    {
        LineT = Lines[n].LTab;
        m++;
    }
    bool *ES = new bool[m];
     for (n=0, m=-1, LineT = 0; n<N; n++)
    {
        if (Lines[n].err == 0.0)
        {
            QMessageBox::information(this, "MolSpektAnalysis", 
                                     "Error: The data contains lines/levels with uncertainty equal to 0!");
            nES = nGS = 0;
            return;
        }
        if (Lines[n].LTab != LineT || (n>0 ? Lines[n-1].isTE != Lines[n].isTE : true))
        {
            LineT = Lines[n].LTab;
            if (Lines[n].isTE)
            {
                nES++;
                ES[++m] = true;
            }
            else 
            {
                if (!OnlyTerm) nGS++;
                ES[++m] = false;
            }
        }
        if (Lines[n].vss <= maxv && Lines[n].Jss <= maxJ) DL[m]++;
    }
    //int mav[maxJ + 1];
    //for (n=0; n <= maxJ; n++) mav[n] = maxv; 
    mJ = new int[NIso];
    mv = new int*[NIso];
    if (nES > 0)
    {
        nEL = new int[nES];
        EL = new TermEnergy*[nES];
    }
    else
    {
        nEL = 0;
        EL = 0;
    }
    if (nGS > 0)
    {
        nGP = new int[nGS];
        nGL = new int[nGS];
        PS = new int*[nGS];
        GL = new TableLine*[nGS];
    }
    else
    {
        nGP = 0;
        nGL = 0;
        PS = 0;
        GL = 0;
    }
    int lPN, PStart[MaxProgressions];
    for (n=0; n < NIso; n++) mJ[n] = -1;
    for (n=e=g=j=0; n<=m; n++)
    {
        if (ES[n]) 
        {
            EL[e] = new TermEnergy[DL[n]];
            for (i=0; i < DL[n]; j++) if (Lines[j].vss <= maxv && Lines[j].Jss <= maxJ)
            {
                EL[e][i].DevR = Lines[j].DevR;
                EL[e][i].E = Lines[j].WN;
                EL[e][i].Iso = Lines[j].Iso;
                EL[e][i].J = Lines[j].Jss;
                EL[e][i].dev = Lines[j].dev;
                EL[e][i].ef = (Lines[j].Jss != Lines[j].Js ? true : false);
                EL[e][i].err = Lines[j].err;
                EL[e][i].v = Lines[j].vss;
                EL[e][i].FC = Lines[j].FC;
                EL[e][i].nDig = Lines[j].nDig;
                EL[e][i].Row = Lines[j].Row;
                EL[e][i].WeightFactor = Lines[j].WeightFact;
                if (EL[e][i].J > mJ[EL[e][i].Iso]) mJ[EL[e][i].Iso] = EL[e][i].J;
                i++;
            }
            nEL[e++] = DL[n];
        }
        else if (!OnlyTerm)
        {
            GL[g] = new TableLine[DL[n]];
            for (i=0, lPN = -1, nGP[g] = 0; i < DL[n]; j++) if (Lines[j].vss <= maxv && Lines[j].Jss <= maxJ)
            {
                GL[g][i] = Lines[j];
                if (lPN != GL[g][i].PN)
                {
                    lPN = GL[g][i].PN;
                    PStart[nGP[g]++] = i;
                }
                if (GL[g][i].Jss > mJ[GL[g][i].Iso]) mJ[GL[g][i].Iso] = GL[g][i].Jss;
                i++;
            }
            PS[g] = new int[nGP[g]];
            for (i=0; i < nGP[g]; i++) PS[g][i] = PStart[i];
            nGL[g++] = DL[n];
        }
    }
    //printf("Ende Schleife\n");
    for (n=0; n < NIso; n++) 
    {
        //printf("NIso=%d, mJ[%d]=%d\n", NIso, n, mJ[n]);
        if (mJ[n] > -1) 
        {
            mv[n] = new int[mJ[n] + 1];
            for (i=0; i <= mJ[n]; i++) mv[n][i] = -1;
        }
        else mv[n] = 0;
    }
    //printf("Vor mv, nES=%d, nGS=%d\n", nES, nGS);
    for (n=Mv=0; n < nES; n++) for (i=0; i < nEL[n]; i++) if (EL[n][i].v > mv[EL[n][i].Iso][EL[n][i].J])
    {
        mv[EL[n][i].Iso][EL[n][i].J] = EL[n][i].v;
        if (EL[n][i].v > Mv) Mv = EL[n][i].v;
    }
    for (n=0; n < nGS; n++) for (i=0; i < nGL[n]; i++) 
            if (GL[n][i].vss > mv[GL[n][i].Iso][GL[n][i].Jss])
    {
        mv[GL[n][i].Iso][GL[n][i].Jss] = GL[n][i].vss;
        if (GL[n][i].vss > Mv) Mv = GL[n][i].vss;
    }
    //printf("Vor Delete\n");
    delete[] Lines;
    //printf("Ende getFitData\n");
    Worker->setFitData(mJ, mv, Mv, nES, nEL, EL, nGS, nGP, nGL, PS, GL, ES, MaxJ, Maxv);
}

void Potential::getCoupledPotentials(int &nCoupled, Potential** &potentials, int* &components)
{
    nCoupled = NCoupled;
    potentials = CoupledPot;
    components = CoupledComp;
}

CoupledSineWaveFunc* Potential::getCoupledSineWaveFunc()
{
    return WaveFuncs;
}

void Potential::setCouplingData(CoupledSineWaveFunc* WF, int nCoupled, Potential** coupledPot, int* coupledComp)
{
    if (NCoupled  > 0)
    {
        delete[] CoupledPot;
        delete[] CoupledComp;
        if (WaveFuncs != 0 ? WaveFuncs->Delete() : false) delete WaveFuncs;
    }
    WaveFuncs = WF;
    WF->Assign();
    NCoupled = nCoupled;
    CoupledPot = coupledPot;
    CoupledComp = coupledComp;
}

void Potential::setCouplingData(CoupledSineWaveFunc* WF, int nCoupled, ElState** CoupledStates, int* coupledComp)
{
    int n, m, *CC;
    Potential **Pot = new Potential*[nCoupled];
    if (WaveFuncs != 0 ? WaveFuncs->Delete() : false) delete WaveFuncs;
    WaveFuncs = WF;
    WF->Assign();
    if (NCoupled != nCoupled - 1)
    {
        if (NCoupled > 0) delete[] CoupledComp;
        CoupledComp = new int[nCoupled];
    }
    for (n=0, m=-1; n < nCoupled; n++)
    {
        if (n>0 ? CoupledStates[n] != CoupledStates[n-1] : false) m++;
        Pot[n] = (m>=0 ? CoupledPot[m] : this);
        CoupledComp[n] = coupledComp[n]; 
    }
    NCoupled = nCoupled;
    delete[] CoupledPot;
    CoupledPot = Pot;
    for (n=0; n < NCoupled; n++)
    {
        Pot = new Potential*[NCoupled];
        CC = new int[NCoupled];
        for (m=0; m < nCoupled; m++)
        {
            Pot[m] = CoupledPot[m];
            CC[m] = coupledComp[m];
        }
        if (CoupledPot[n] != this) CoupledPot[n]->setCouplingData(WF, NCoupled, Pot, CC);
        if (!CoupledPot[n]->isAssigned()) CoupledStates[n]->addPotential(CoupledPot[n]);
    }
}

double Potential::getD0(int NumWFPoints)
{
    int I = Worker->getRefIso(), nv = 1;
    double h = (rmax - rmin) / (NumWFPoints - 1);
    double *U = getPoints(rmin, rmax, NumWFPoints), *IsoF = Worker->getIsoF_ForReading(), MF = IsoF[I] * h*h;
    double **WF = Create(1, NumWFPoints), *Y = new double[NumWFPoints];
    printf("Nach Create WF\n");
    //aPot->GetMin(R, M);
    double F, SQ, E;
    int CP, ri, ra;
    NumerovCooley(Y, MF, nv, &E, WF, &ri, &ra, &CP, &F, &SQ, NumWFPoints, U);
    delete[] U;
    delete[] Y;
    Destroy(WF, 1);
    return Worker->getUinf() - E;
}

bool Potential::getCalcDiagFuncs()
{
    return Worker->getCalcDiagFuncs();
}

void Potential::setCalcDiagFuncs(bool calc)
{
    if (calc && !Worker->getCalcDiagFuncs()) 
    {
        Worker->setCalcDiagFuncs(calc);
        getFQS(true);
    }
    else Worker->setCalcDiagFuncs(calc);
}

void Potential::getDiagFuncs(double *&WaveFSum, double *&rWWFS, double &RStart, double &h, int &NP, int NumWFPoints)
{
    //printf("Potential::getDiagFuncs\n");
    Worker->getDiagFuncs(WaveFSum, rWWFS, NumWFPoints);
    RStart = rmin;
    h = (rmax - rmin) / (NumWFPoints - 1);
    NP = NumWFPoints;
}

void Potential::showBadList(int NumWFPoints)
{
    int NR, NC;
    QString **Data;
    if (!Fit->isRunning() && !Worker->isFitResultAv()) setFitData();
    Worker->showBadList(Data, NR, NC, NumWFPoints);
    if (NR == 0 && !Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
        "The badlist is empty, because all lines are explained sufficiently well by the current potential.");
        Worker->unlockBadList();
        return;
    }
    if (BadList == 0) 
    {
        BadList = MW->createTableWindow();
        BadList->setType(TextTable1);
        BadList->setWindowTitle("Badlist for the potential " + getName());
        connect(Worker, SIGNAL(badListChanged()), this, SLOT(badListChanged()));
    }
    if (NC == 7) BadList->setHorizontalHeader(QStringList() << "Iso" << "v" << "J" << "Energy" 
                    << "Uncertainty" << "Obs-calc" << "Dev ratio");
    else BadList->setHorizontalHeader(QStringList() << "Iso" << "v''" << "J''" << "v'" << "J'" 
                << "Energy" << "Uncertainty" << "Obs-calc" << "Dev ratio");
    BadList->setData(Data, NR, NC);
    BadList->setEditable(false);
    Worker->unlockBadList();
    if (!BadList->isVisible()) BadList->show();
    connect(BadList, SIGNAL(TabRowChanged(QString*,int)), this, SLOT(badListChanged(QString*,int)));
}

void Potential::startFit()
{
    if (Worker->isSplinePot()) FitSplinePot();
    else if (Worker->isAnaPotAv()) FitAnaPot();
}

void Potential::badListChanged(double aFQS, double aLambda, double mLSFQS)
{
    if (FQSv == -1.0) 
    {
        bFQS = FQSv = aFQS;
        NFitIt = 0;
    }
    else if (bFQS > aFQS) bFQS = aFQS;
    if (LogWindow != 0)
    {
        QString Text = QString("%1 %4: FQS=%2, bFQS=%3").arg(NFitIt++).arg(aFQS, 0, 'g', 12).arg(bFQS, 0, 'g', 12);
        LogWindow->addTextRow(aLambda > 0 ? Text.arg("LM").append(", l=%1, LS=%2").arg(aLambda, 0, 'g', 1).arg(mLSFQS, 0, 'g', 3) 
                                          : Text.arg("SVD"));
    }
}

void Potential::badListChanged(QString* Items, int N)
{
    if (Fit->isRunning()) return;
    TableLine L;
    L.Iso = Items[0].toInt();
    L.Jss = Items[2].toInt();
    L.vss = Items[1].toInt();
    if (N==7)
    {
        L.dev = Items[5].toDouble();
        L.DevR = Items[6].toDouble();
        L.err = Items[4].toDouble();
        L.WN = Items[3].toDouble();
    }
    else
    {
        L.dev = Items[7].toDouble();
        L.DevR = Items[8].toDouble();
        L.err = Items[6].toDouble();
        L.WN = Items[5].toDouble();
        if (!Items[3].isEmpty()) 
        {
            L.Js = Items[4].toInt();
            L.vs = Items[3].toInt();
        }
    }
    fitData->updateRow(&L);
    Worker->badListChanged(Items, N, L);
}

double Potential::getFQS(int NumWFPoints, bool calc, double* StdDev, double *sigma, double *FQS_PAL)
{
    if (molecule == 0) return -1.0;
    if (!Fit->isRunning()) setFitData();
    double rFQS = Worker->getFQS(NumWFPoints, calc, StdDev, false, 0.0, sigma, FQS_PAL), *dev, *RWErr;
    int N, *Row;
    Worker->getDeviations(N, Row, dev, RWErr);
    fitData->setDev(dev, Row, N);
    Worker->unlockBadList();
    return rFQS;
}

void Potential::calcFQS_SFQS(bool ****SFQSU, double SFQSRad)
{
    if (molecule == 0 || Fit->isRunning()) 
    {
        printf("calcFQS_SFQS: invalid conditions for starting!\n");
        return;
    }
    if (SFQSU != 0) Worker->setSFQSU(SFQSU);
    setFitData();
    showFitResult = false;
    Fit->CalcFQS(SFQSRad);
}

double Potential::CreateSplineFromAnalyticalPot()
{
    Potential *NPot = MW->CreatePotential();
    if (NPot == 0) return 0.0;
    double MinR, MinE, MaxE, E, R, RMin, RMax, *LRC, A, iO, iA, iExp;
    SplinePoint *R1, *R2;
    int n, N=0, N1, N2 = 0, *PLRC, NLRC;
    QString F;
    QDialog *D = new QDialog(this);
    QGridLayout *L = new QGridLayout(D);
    QLineEdit *NPE = new QLineEdit("20", D);
    QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
    D->setWindowTitle("Create spline from analytical potential");
    L->addWidget(new QLabel("Number of points:", D), 0, 0);
    L->addWidget(NPE, 0, 1);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(OK, 2, 0);
    L->addWidget(Cancel, 2, 1);
    connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
    while (N < 2 || N > 100)
    {
        if (D->exec() == QDialog::Rejected)
        {
            delete D;
            return 0.0;
        }
        N = NPE->text().toInt();
        if (N < 2 || N > 100) QMessageBox::information(this, "MolSpektAnalysis", "You must choose a valid integer between 2 and 100");
    }
    delete D;
    AnaPot *aPot = Worker->getAnaPotForReading();
    Worker->getLRCoeffForReading(NLRC, PLRC, LRC, &RMax, true);
    aPot->getInnerWall(RMin, iO, iExp, iA);
    Worker->unlock();
    aPot->getMinimum(MinR, MinE);
    MaxE = getPoint(RMax);
    if ((E = Worker->getPoint(RMin)) > MaxE) MaxE = E;
    R1 = new SplinePoint[N];
    R2 = new SplinePoint[N];
    R = 0.30 * (RMax - RMin) / (N-1);
    R2[0].x = MinR + R;
    R2[0].y = Worker->getPoint(R2[0].x);
    for (n=0, E = (MaxE - MinE); n<N-1; 
            R2[n].y = Worker->getPoint(R2[n].x = R2[n-1].x + 1.1 * E * R / (0.1 * E + MaxE - R2[n-1].y))) 
    {
        if (R2[n].x < RMax) N2 = n + 2;
        n++;
    }
    E = (MinR - RMin) / R + 0.5;
    N1 = int(ceil(E * double(N) / (E + double(N2))));
    for (R1[n=0].x = MinR - 0.5 * (R = (MinR - RMin) / (double(N1) - 0.5)); n < N1; n++, R1[n].x = R1[n-1].x - R) 
        R1[n].y = Worker->getPoint(R1[n].x);
    for (R = fabs((R1[0].y - (E = R2[0].y)) / (R2[0].y - MinE) * ((A = R2[0].x) - MinR)); 
            fabs(R1[0].y - R2[0].y) > 0.01 * (R1[0].y - MinE); 
            R2[0].y = Worker->getPoint(R2[0].x += (R1[0].y > (E = R2[0].y) ? R : -R))) 
         if ((E - R1[0].y) * (R2[0].y - R1[0].y) < 0.0) R/=2; 
    for (n=1, N2 = N - N1, A -= R2[0].x; n < N2; n++) R2[n].y = Worker->getPoint(R2[n].x -= A);
    R = ((R2[N2 - 1].x - RMax) * R2[0].y) / (RMax * (R2[N2 - 1].y - R2[0].y));
    for (n=1; n < N2; n++) R2[n].y = Worker->getPoint(R2[n].x *= (R2[0].y / (R2[0].y + R * (R2[n].y - R2[0].y))));
    SplinePoint *points = new SplinePoint[N];
    for (n=0; n < N1; n++) points[n] = R1[N1 - n - 1];
    for ( ; n<N; n++) points[n] = R2[n - N1];
    delete[] R1;
    delete[] R2;
    NPot->setName(getName() + "spline");
    NPot->setSource("Analytical potential from " + getSource());
    F = getFileName();
    NPot->setFileName(F.right(F.length() - F.indexOf('.') - 1) + "spline.pot");
    NPot->setWorker(new SplinePot(NPot->Fit));
    NPot->Worker->setUinf(aPot->getUinf());
    int *NPLRC = new int[NLRC];
    double *nLRC = new double[NLRC];
    for (n=0, R=0.0; n < NLRC; n++) 
    {
        NPLRC[n] = PLRC[n];
        nLRC[n] = LRC[n];
        R += LRC[n] * double(PLRC[n]) * pow(points[N-1].x, -PLRC[n] - 1);
    }
    if (Worker->isAdCorrA())
    {
        double *adCorr, RIso1, RIso2, b, Rm;
        int NAdCorr, PAdCorr, TAdCorr;
        bool *AdCorrFree;
        Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, &AdCorrFree);
        NPot->Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, AdCorrFree);
    }
    NPot->Tab->setHorizontalHeaderLabels(QStringList() << "R [A]" << "E [cm^-1]" << "d^2E/dR^2");
    NPot->Tab->setColumnCount(3);
    NPot->Worker->setSplinePot(N, points, NLRC, NPLRC, nLRC, iA, iO, iExp, R);
    NPot->UpdateTab();
    State->addPotential(NPot);
    NPot->show();
    //return NPot->FitSplinePot(true, this);
    return 0.0;
}

void Potential::FitSplinePot(bool showResult, int MaxIt, double threshFakt)
{
    if (!Worker->isSplinePot()) 
    {
        if (Worker->isAnaPotAv())
        {
            CreateSplineFromAnalyticalPot();
            return;
        }
        QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                 tr("Error: no spline potential available!"));
        return;
    }
    if (molecule == 0) 
    {
        QMessageBox::information(this, "QT4MolSpektAn", 
                                 "Error: potential is not assigned to a molecule!");
        return;
    }
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
                        "A fit of the potential cannot be started while there is already a fit of the same potential running!");
        return;
    }
    mWasMoving = false;
    setFitData();
    if (!Worker->isFitDataAv()) return; 
    showFitResult = showResult;
    FQSv = -1.0;
    if (LogWindow == 0) initLogWindow();
    LogWindow->setWindowTitle("Fit of spline potential " + getName());
    Fit->FitSplinePot(MaxIt, threshFakt);
}

void Potential::MonteCarloSimIteration(double *MCSTab, double UncFact)
{
    if (Fit->isRunning()) return;
    setFitData();
    if (LogWindow == 0) initLogWindow();
    LogWindow->setWindowTitle("MCS Fit Potential " + getName());
    showFitResult = false;
    FQSv = -1.0;
    Fit->RunMonteCarloIteration(MCSTab, UncFact);
}

void Potential::CreateFirstPotFitTraceTabRows(FitResult Result, double FQS)
{
    QStringList Row;
    int c, NC = Tab->rowCount() - 2;
    for (c=0; c < NC; c++) Row << Tab->item(c, 0)->text();
    Row << "FQS" << "FQS_Bad" << "N_Bad" << "N_Bad_PAL" << "Sigma";
    PotFitTraceTab.clear();
    PotFitTraceTab << Row.join("\t");
    for (c=0, Row.clear(); c < NC; c++) Row << Tab->item(c, 1)->text();
    Row << QString::number(FQS, 'g', 10) << QString::number(Result.FQS_Bad, 'g', 10) << QString::number(Result.NBad) << QString::number(Result.NBadPAL) 
        << QString::number(Result.Sigma, 'f', 10);
    PotFitTraceTab << Row.join("\t");
}

void Potential::potentialImproved(PotentialData *Data)
{
    QString FName = getFileName();
    QString TFile = FName.left(FName.lastIndexOf('.') + 1) + "temp.pot";
    if (fitData != 0)
    {
        int N, *Row;
        double *dev, *RWErr;
        Worker->getDeviations(N, Row, dev, RWErr);
        fitData->setDev(dev, Row, N);
        if (RWErr != 0) fitData->setRWErr(RWErr);
        Worker->unlockBadList();
    }
    if (BadList != 0)
    {
        int NR, NC;
        QString **Data;
        Worker->getBadList(Data, NR, NC);
        if (BadList->getNumColumns() != NC)
        {
            if (NC == 7) BadList->setHorizontalHeader(QStringList() << "Iso" << "v" << "J" << "Energy" 
                                                        << "Uncertainty" << "Obs-calc" << "Dev ratio");
            else BadList->setHorizontalHeader(QStringList() << "Iso" << "v''" << "J''" << "v'" << "J'" 
                                                << "Energy" << "Uncertainty" << "Obs-calc" << "Dev ratio");
        }
        BadList->setData(Data, NR, NC);
        BadList->setEditable(false);
        Worker->unlockBadList();
    }
    FitResult Result = Data->fitResult;
    double FQS = Data->FQS;
    UpdateTab(Data);
    if (writePotFitTraceTab)
    {
        QStringList Row;
        int c, NC = Tab->rowCount() - 2;
        for (c=0, Row.clear(); c < NC; c++) Row << Tab->item(c, 1)->text();
        Row << QString::number(FQS, 'g', 10) << QString::number(Result.FQS_Bad, 'g', 10) << QString::number(Result.NBad) << QString::number(Result.NBadPAL) 
            << QString::number(Result.Sigma, 'f', 10);
        PotFitTraceTab << Row.join("\t");
    }
    writeData(TFile, false);
    setFileName(FName);
    Changed();
}

void Potential::fitFinished()
{
    double FQS, WeightSum, *WeightF;
    int NWeightF, *WeightFRow, NumFreePar;
    Worker->getFitResult(FQS, WeightSum, NWeightF, WeightFRow, WeightF, NumFreePar);
    if (NWeightF > 0) fitData->setWeightFactors(NWeightF, WeightFRow, WeightF);
    Worker->unlockBadList();
    if (Fit->getTask() == PotFit::calcFQS)
    {
        emit SFQScalculated(threadNum, Worker->getSFQS(), FQS);
        return;
    }
    if (!Worker->isLevenbergMarquardtRunning()) potentialImproved(Worker->getPotentialData());
    if (showFitResult) QMessageBox::information(this, tr("Fit result"), 
                                 "Initial chiSq=" + QString::number(FQSv, 'g', 6) + ", final chiSq="
                                 + QString::number(FQS, 'g', 6) + ", Sigma=" 
                                 + QString::number(sqrt(FQS/(WeightSum - double(NumFreePar)))));
    UpdateTab();
    Changed();
    if (Fit->getTask() == PotFit::MonteCarloSimulation) 
    {
        int N, *Rows;
        double *E = Worker->getMCSE(N, Rows);
        fitData->setMCSEnergies(E, Rows, N);
    }
    if (getTexTableSplineFitRunning) getTexTable(FQS);
    emit FitFinished(threadNum, FQS);
}

void Potential::fitTerminated()
{
    emit FitTerminated(threadNum);
}

void Potential::restartFit()
{
    if (Worker->isFitStopped())
    {
        if (Fit->isRunning()) Fit->wait();
        Fit->restartFit();
    }
    else if (Fit->isRunning()) Worker->stopFit(false);
    Worker->unlock();
}

void Potential::cancelFit(bool wait)
{
    Worker->stopFit(true);
    Fit->terminate();
    if (wait && Fit->isRunning()) Fit->wait();
}

void Potential::initLogWindow()
{
    LogWindow = new PotFitLogWindow();
    connect(LogWindow, SIGNAL(CancelFit(bool)),this, SLOT(cancelFit(bool)));
    connect(LogWindow, SIGNAL(StopFit()), this, SLOT(stopFit()));
    connect(LogWindow, SIGNAL(StartFit()), this, SLOT(startFit()));
    connect(LogWindow, SIGNAL(RestartFit()), this, SLOT(restartFit()));
    connect(Worker, SIGNAL(fitRunning()), LogWindow, SLOT(FitRunning()));
    connect(Worker, SIGNAL(fitStopped()), LogWindow, SLOT(FitStopped()));
    connect(this, SIGNAL(FitFinished(int, double)), LogWindow, SLOT(FitFinished()));
    connect(this, SIGNAL(FitTerminated(int)), LogWindow, SLOT(FitTerminated()));
    MW->showMDIChild(LogWindow);
    MW->setAskForQuit(LogWindow);
}

void Potential::FitAnaPot(bool showResult, bool robustWeighting, bool UseSVD, double threshhold, bool transform, bool ask,
                          bool useLeveMarq, int MaxIt, bool ATEWF, double hCi_cP, bool WritePotFitTraceTab)
{
    int NC, numPoints = 10000, nLRC;
    double Ri, Ra;
    double Rmin, Rmax, iEx = Worker->getIExp(), bmin = -0.1, bmax = -0.59, bstep = -0.01, A_ex = 0.0, beta, gamma;
    bool improvePotential = ((Worker->isAnaPotAv() || Worker->isSplinePot()) && !transform);
    int *pLRC = Worker->getSplinePLRC();
    PotentialType FType;
    mWasMoving = false;
    if (!WritePotFitTraceTab && writePotFitTraceTab) disconnect(Worker, SIGNAL(firstFCFCalculated(FitResult,double)), this, SLOT(CreateFirstPotFitTraceTabRows(FitResult,double)));
    else if (!writePotFitTraceTab && WritePotFitTraceTab) connect(Worker, SIGNAL(firstFCFCalculated(FitResult,double)), this, SLOT(CreateFirstPotFitTraceTabRows(FitResult,double)));
    writePotFitTraceTab = WritePotFitTraceTab;
    Worker->getSplinePotPointandCoeffNumbers(NC, nLRC);
    Worker->getSplineMinMaxR(Ri, Ra);
    if (Ri == 0.0)
    {
        Ri = rmin + 1.0;
        Ra = 10.0;
    }
    if (ask)
    {
        if (Worker->isAnaPotAv())
        {
            double *aLRC, aiC, aiO;
            int N, n;
            AnaPot *aPot = Worker->getAnaPotForReading();
            bool *LRCFree = aPot->getLRCFree();
            Worker->getLRCoeffForReading(nLRC, pLRC, aLRC, &Ra, true);
            aPot->getInnerWall(Ri, aiO, iEx, aiC);
            aPot->getNC(N, NC, nLRC);
            if (hCi_cP <= 0.0) for (hCi_cP = -1.0, n=0; n < nLRC; n++) if (LRCFree[n] && pLRC[n] >= 12) hCi_cP = 0.0;
            Worker->unlock();
        }
        else NC = 0;
        FitAnaPotDialog *D = new FitAnaPotDialog(improvePotential, Ri, Ra, iEx, NC, nLRC, pLRC, hCi_cP, Worker->getType(), this);
        if (!D->exec(improvePotential, useLeveMarq, UseSVD, Ri, Ra, iEx, NC, nLRC, pLRC, Rmin, Rmax, numPoints, 
            bmin, bmax, bstep, ATEWF, hCi_cP, FType)) 
        {
            delete D;
            return;
        }
        delete D;
        if (bstep == 0.0) bstep = -0.01;
        if (bstep * (bmax - bmin) < 0.0)
        {
            double db = bmax;
            bmax = bmin;
            bmin = db;
        }
    }
    else
    {
        Rmin = Ri;
        Rmax = Ra;
        NC = 15;
    }
    if (!improvePotential)
    {
        int n;
        double *Sig, *U;
        QString F;
        Potential *nPot = (!transform ? MW->CreatePotential() : this);
        AnaPot *TAPot = new AnaPot(nPot->Fit);
        double h = (Rmax - Rmin) / (numPoints - 1), b, Rm, bFQS = 1e50, FQS;
        double uinf = Worker->getUinf();
        U = Worker->getPoints(Rmin, Rmax, numPoints);
        Sig = new double[numPoints];
        if (Worker->isAnaPotAv()) 
        {
            AnaPot *aPot = Worker->getAnaPotForReading();
            aPot->getExchangeCoeff(A_ex, beta, gamma);
            if (A_ex != 0.0) TAPot->setExchangeCoeff(A_ex, beta, gamma, aPot->isExcIntFree());
            Worker->unlock();
        }
        int *PLRC;
        double *LRC;
        bool *LRCFree = Worker->getLRCFree(), *lrcFree = new bool[nLRC];
        for (n=0; n < nLRC; n++) lrcFree[n] = LRCFree[n];
        Worker->getLRCoeffForReading(nLRC, PLRC, LRC);
        TAPot->setLRCoeff(Ra, nLRC, PLRC, LRC, lrcFree);
        for (n=0; n < numPoints; n++) Sig[n] = 1e-4;
        getMinimum(Rm, FQS);
        for (b = bmin; (bstep < 0 ? b > bmax + bstep : b < bmax + bstep); b+= bstep) 
            if ((FQS = TAPot->fitPotential(numPoints, Rmin, h, U, Sig, NC, nLRC, pLRC, Ri, Ra, Rm, b, uinf, iEx)) < bFQS)
        {
            bFQS = FQS;
            nPot->setWorker(new AnaPot(*TAPot));
        }
        if (Worker->isAdCorrA() && !transform)
        {
            int NAdCorr, TAdCorr, PAdCorr;
            double *adCorr, Rm, b, RIso1, RIso2;
            bool *adCorrFree;
            Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, &adCorrFree);
            nPot->Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, adCorrFree);
        }
        delete[] U;
        delete[] Sig;
        nPot->UpdateTab();
        nPot->setName(getName() + "AnaPot");
        nPot->setSource("Spline potential from " + getSource());
        F = getFileName();
        nPot->setFileName(F.right(F.length() - F.indexOf('.') - 1) + "AnaPot.pot");
        if (!transform)
        {
            if (State != 0) State->addPotential(nPot);
            nPot->show();
        }
        //if (!useLeveMarq && !UseSVD) return bFQS;
        if (!transform) nPot->FitAnaPot(showResult, robustWeighting, UseSVD, threshhold, false, false, useLeveMarq, MaxIt);
        return;
    }
    else
    {
        if (Worker->getNIso() == 0)
        {
            QMessageBox::information(this, "MolSpektAnalysis", "The potential has to be assigned to a molecule with atomic data first!");
            return;
        }
        if (Fit->isRunning())
        {
            QMessageBox::information(this, "MolSpektAnalysis", 
                            "A fit of the potential cannot be started while there is already a fit of the same potential running!");
            return;
        }
        if (LogWindow == 0) initLogWindow();
        LogWindow->setWindowTitle("Fit of analytical potential " + getName());
        setFitData(0, 0, -1);
        showFitResult = showResult;
        FQSv = -1.0;
        Fit->FitAnaPot(robustWeighting, UseSVD, threshhold, useLeveMarq, MaxIt, ATEWF, hCi_cP);
        Fit->start();
    }
    
    /*if (State != 0 ? (fitData = State->getFitData()) == 0 : true)
    {
        QMessageBox::information(this, tr("MolSpektAnalysis"), 
                             "Error: No fit data available!");
        return;
    }
    if (aPot == 0 && LPot == 0)
    {
        QMessageBox::information(this, tr("MolSpektAnalysis"), 
                             "Error: No initial potential data available!");
        return;
    }
    int NC, nLRC, NL, i, o, n, *lp, *pLRC, p, q, NS;
    double *C, *lRC, Rm, Ra, b, Tm, De;
    TableLine *data, *data1, dB;
    LineTable *LSource = 0, *Sources[MaxLineTables];
    fitData->getData(data1, NL);
    for (n = NS = i = 0; n < NL; n++) 
    {
        if (!data1[n].isTE)
        {
            if (data1[n].LTab != LSource)
            {
                LSource = data1[n].LTab;
                for (i=0; (i < NS ? LSource != Sources[i] : false); i++) ;
                if (i == NS) Sources[NS++] = LSource;
            }
            data1[n].SourceN = i;
        }
        else data1[n].SourceN = -1;
    }
    for (i=1; i < NL; i++) 
        for (o=((n=i)-1)/2; 
             (o>=0 ? data1[n].Iso < data1[o].Iso || (data1[n].Iso == data1[o].Iso && data1[n].Jss < data1[o].Jss) 
                   : false); o=((n=o)-1)/2)
    {
        dB = data1[n];
        data1[n] = data1[o];
        data1[o] = dB;
    }
    data = new TableLine[NL];
    for (i=0; i < NL; i++)
    {
        data[i] = data1[0];
        for (o=1, n=0; o < NL; o=2*(n=o)+1)
        {
            if (o+1 < NL)
                if (data1[o].Iso > data1[o+1].Iso || (data1[o].Iso == data1[o+1].Iso && data1[o].Jss > data1[o+1].Jss)) 
                    o++;
            if (data1[o].Iso >= 1000) break;
            data1[n] = data1[o];
            data1[o].Iso += 1000;
        }
    }
    delete[] data1;
    if (aPot != 0) 
    {
        aPot->getPotCoeff(b, Rm, Tm, De, NC, C);
        aPot->getLRCoeff(Ra, nLRC, lp, lRC);
    }
    else 
    {
        LPot->getNC(i, NC, nLRC);
        NC += 3;
    }
    MrqLev ML(NC + nLRC, NL);
    if (aPot != 0) 
    {
        for (i=0; i < NC; i++)
        {
            ML.a[i] = C[i];
            ML.ia[i] = true;
        }
        for (n=0; n < nLRC; n++)
        {
            ML.a[i] = lRC[n];
            ML.ia[i++] = true;
        }
    }
    else 
    {
        LPot->getCoefficients(p, q, ML.a[0], ML.a[1], ML.a[2], Rm, NC, C, nLRC, lp, lRC);
        ML.ia[0] = ML.ia[1] = ML.ia[2] = true;
        for (n=0, i=3; n < NC; n++)
        {
            ML.a[i] = C[n];
            ML.ia[i++] = true;
        }
        for (n=0; n < nLRC; n++)
        {
            ML.a[i] = lRC[n];
            ML.ia[i++] = false;
        }
    }
    delete[] lRC;
    delete[] C;
    calcDerivatives(data, NL, ML.E, ML.dyda);
    for (n=0, i=-1; n < NL; n++) 
    {
        ML.sig[n] = data[n].err;
        ML.y[n] = data[n].WN;
        if (!data[n].isTE && data[n].PN > i) i = data[n].PN;
    }
    if (i>=0)
    {
        double **SE = Create(NS, i+1), **WS = Create(NS, i+1), W;
        for (o=0; o < NS; o++) for (n=0; n<=i; n++) SE[o][n] = WS[o][n] = 0.0;
        for (n=0; n < NL; n++) if (!data[n].isTE)
        {
            W = 1.0 / (data[n].err * data[n].err);
            WS[data[n].SourceN][data[n].PN] += W;
            SE[data[n].SourceN][data[n].PN] += W * data[n].WN;
        }
        for (o=0; o < NS; o++) for (n=0; n<=i; n++) if (WS[o][n] != 0.0) SE[o][n] /= WS[o][n];
        for (n=0; n < NL; n++) if (!data[n].isTE) ML.y[n] -= SE[data[n].SourceN][data[n].PN];
        Destroy(SE, NS);
        Destroy(WS, NS);FileName.left(FileName.length() - 4) + "_i.pot"
    }
    ML.init();
    printf("init: %e ma=%d, mfit=%d\n", ML.chisq, ML.ma, ML.mfit);
    for (i=0; i < 1000; i++)
    {
        pLRC = new int[nLRC];
        lRC = new double[nLRC];
        C = new double[NC];
        if (LPot != 0)
        {
            for (n=0, o=3; n < NC; n++) C[n] = ML.atry[o++];
            for (n=0; n < nLRC; n++)
            {
                pLRC[n] = lp[n];
                lRC[n] = ML.atry[o++];
            }
            LPot->setCoefficients(p, q, ML.atry[0], ML.atry[1], ML.atry[2], Rm, NC, C, nLRC, pLRC, lRC);
        }
        else if (aPot != 0)
        {
            for (o=0; o < NC; o++) C[o] = ML.atry[o];
            for (n=0; n < nLRC; n++)
            {
                pLRC[n] = lp[n];
                lRC[n] = ML.atry[o++];
            }
            aPot->setPotCoeff(b, Rm, Tm, De, NC, C);
            aPot->setLRCoeff(Ra, nLRC, pLRC, lRC);
        }
        calcDerivatives(data, NL, ML.E, ML.dyda);
        if (!ML.fit()) 
        {
            printf("finished: chisq=%e\n", ML.chisq);
            break;
        }
        if (ML.alamda > 1e300)
        {
            printf("lamda to large!\n");
            break;
        }
        printf("%d: %e nd=%d, l=%e\n", i, ML.chisq, ML.done, ML.alamda);
    }
    delete[] data;
    UpdateTab();
    Changed();*/
    
    /*int n, NP = NumWFPoints, NK = 5, np, ri, pi, pa, nP, NMax = 0;
    double RMin, RMax, Max, *U, Tm, b, Rm, bb, Err, sErr, De, rMin = rmin, rMax = rmax;
    if (aPot != 0) 
    {
        aPot->getConR(RMin, RMax);
        De = aPot->GetDe();
    }
    else if (points != 0)
    {
        RMin = points[0].x;
        RMax = points[numSplinePoints - 1].x;
        De = points[numSplinePoints - 1].y;
    }
    else
    {
        QMessageBox::information(this, tr("MolSpektAnalysis"), 
                             "Error: No initial potential data available!");
        return;
    }
    double h = (rMax - rMin) / double(NP - 1), r;
    setBlockChangeSignal(true);
    for (r = rMin, pi=0; r <= RMin; r+=h, pi++) ;
    for (RMin = r, pa = pi; r < RMax; r+=h, pa++) ;
    RMax = r-h;
    nP = (pa--) - pi;
    double *R = new double[nP], M, MF = 0.0833333333333 * h * h;
    U = getPoints(rMin, rMax, NP);
    for (R[0] = RMin, n=1; n < nP; n++) R[n] = R[n-1] + h;
    for (Max = U[NP - 1], n = pi; n < pa; n++) if (U[n] > Max && U[n] > U[n-1] && U[n] > U[n+1]) 
            Max = U[n];
    for (n = pi; (n < pa ? U[n] > Max : false); n++) ;
    if (n == pa)
    {
        QMessageBox::information(this, tr("MolSpektAnalysis"), 
                            "Error: This fit function should not be used for repulsive potentials!");
        delete[] R;
        delete[] U;
        setBlockChangeSignal(false);
        return;
    }
    ri = n - pi;
    pi = n;
    int *lExp = new int[2];
    double *PotK = new double[NK], *lK = new double[2], *Sig = new double[NP], *WF = new double[NP];
    double *UR = new double[NP];
    int i, j, *mJ, **mv, Mv, nEs, *NEL, nGS, *nGP, *nGL, **PS, I, J, v, ***NL = new int**[NIso];
    TermEnergy **EL;
    TableLine **GL;
    TermTable *TT = State->getTermTable();
    double ****Data = TT->getData(), Ri, Ra, ISq, F, CF;
    int Nv = TT->getMaxv(), NJ = TT->getMaxJ(), Wi, Wc, Wa;
    getFitData(mJ, mv, Mv, nEs, NEL, EL, nGS, nGP, nGL, PS, GL);
    for (I=0; I < NIso; I++) if (mJ[I] != -1)
    {
        NL[I] = new int*[mJ[I] + 1];
        for (J=0; J <= mJ[I]; J++) if (mv[I][J] != -1)
        {
            NL[I][J] = new int[mv[I][J] + 1];
            for (v=0; v <= mv[I][J]; v++) NL[I][J][v] = 0;
        }
    }
    for (i=0; i < nEs; i++) for (j=0; j < NEL[i]; j++) NL[EL[i][j].Iso][EL[i][j].J][EL[i][j].v]++;
    for (i=0; i < nGS; i++) for (j=0; j < nGL[i]; j++) NL[GL[i][j].Iso][GL[i][j].Jss][GL[i][j].vss]++;
    printf("nP=%d, pa=%d\n", nP, pa);
    for (n=pi; n < pa; n++) Sig[n] = 0.0;
    for (I=0; I < NIso; I++) if (mJ[I] != -1) for (M = MF * IsoF[I], J=0; J <= mJ[I]; J++) 
                if (mv[I][J] != -1 && J <= NJ)
    {
        CRotPot(J, I, U, UR, rMin, h, NP);
        printf("I=%d, J=%d, Wi=%d, Wc=%d, Wa=%d\n", I, J, Wi, Wc, Wa);
        for (v=0; v <= mv[I][J]; v++) if (NL[I][J][v] > 0 && v <= Nv && Data[0][I][v][J] != 0.0)
        {
            fastWaveFunc(UR, Data[0][I][v][J], M, h, NP, WF, Ri, Ra, ISq, CF, Wi, Wc, Wa);
            ISq = 1.0 / ISq;
            CF *= CF;
            if (Wa > nP) printf("v=%d, nP=%d, Wa=%d\n", v, nP, Wa);
            for (n = Wi; n < Wa; n++)
            {
                F = WF[n] * WF[n] * (n <= Wc ? ISq : ISq * CF);
                if (F > 1.0) printf("I=%d, J=%d, v=%d, n=%d, F=%g\n", I, J, v, n, F);
                Sig[n] += F * F * double(NL[I][J][v]);
            }
        }
    }
    getMinimum(Rm, Tm);
    printf("pi=%d, pa=%d, nP=%d\n", pi, pa, nP);
    for (n = pi, M = 0.0; n < pa; n++) if (Sig[n] > M) M = Sig[n];// i=n;}
    for (NMax = pa - 1, F = 0.01 * M; Sig[NMax] < F; NMax--) ;
    printf("M=%g\n", M);
    for (n = pi, F = sqrt(1.0 / (M *= 1e-3)); n < pa; n++) 
    {
        if (Sig[n] < M) Sig[n] = F;
        else Sig[n] = sqrt(1.0 / Sig[n]);
    }
    np = pa - pi;
    printf("pi=%d, np=%d\n", pi, np);
    printf("R[ri]=%g, rMin=%g, rMin+pi*h=%g\n", R[ri], rMin, rMin + pi * h);
    for (b=-0.1, sErr = 1e300; b>=-0.6; b-=0.01) 
    {
        if ((Err = fitPotCoeff(np, R + ri, U + pi, NK, Rm, b, PotK, Tm, Sig + pi)) < sErr)
        {
            bb = b;
            sErr = Err;
        }
        printf("b=%e, Err=%g\n", b, Err);
    }
    fitPotCoeff(np, R + ri, U + pi, NK, Rm, bb, PotK, Tm, Sig + pi);
    if (aPot != 0) delete aPot;
    else Tab->setHorizontalHeaderLabels(QStringList() 
                    << "coefficient" << "value" << "error" << "sig. digits");
    aPot = new AnaPot;
    aPot->setPotCoeff(b, Rm, Tm, De, NK, PotK);
    aPot->setLRCoeff(RMax, 0, lExp, lK);
    aPot->sdConnectLR(6);
    aPot->sdConnectIW(6, RMin + double(ri) * h);
    delete[] U;
    U = getPoints(rMin, rMax, NP);
    for (n = NMax, M = U[NP - 1]; (n < NP ? U[n] < M && U[n] > U[n-1] : false); n++) ;
    if ((F = rMin + double(n) * h) < RMax)
    {
        aPot->setLRCoeff(F, 0, lExp, lK);
        aPot->sdConnectLR(6);
    }
    UpdateTab();
    for (n=0; n < numSplinePoints; n++) points[n].y = points[n].yss = 0.0;
    printf("Vor delete\n");
    for (I=0; I < NIso; I++) if (mJ[I] != -1)
    {
        for (J=0; J <= mJ[I]; J++) if (mv[I][J] != -1) delete[] NL[I][J];
        delete[] NL[I];
        delete[] mv[I];
    }
    printf("vor delete NL\n");
    delete[] NL;
    delete[] mJ;
    delete[] mv;
    printf("Vor delete PS\n");
    for (n=0; n < nGS; n++) delete[] PS[n];
    printf("Vor NEL\n");
    delete[] NEL;
    delete[] nGL;
    delete[] nGP;
    delete[] PS;
    printf("Vor delete sig\n");
    delete[] Sig;
    printf("WF\n");
    delete[] WF;
    printf("U\n");
    delete[] U;
    printf("R\n");
    delete[] R;
    printf("UR\n");
    delete[] UR;
    setBlockChangeSignal(false);
    Changed();
    printf("Ende  FitAnapot\n");
    */
}

double Potential::FitTangToenniesPot(bool showResult)
{
    if (!Worker->isTTPot())
    {
        int n, NA, NC, numPoints = 10000, CNC, *pCLRC;
        double Re, Te, *CLRC, b;
        bool ExcInt;
        QString F;
        Potential *nPot = MW->CreatePotential();
        TangToenniesPot *TTPot = new TangToenniesPot(nPot->Fit);
        QDialog *D = new QDialog(this);
        QGridLayout *L = new QGridLayout(D);
        QLineEdit *NPE = new QLineEdit("6", D), *NAE = new QLineEdit("7", D);
        QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
        QCheckBox *ExcIntB = new QCheckBox("Use exchange interaction", D);
        D->setWindowTitle("Create Tang-Toennies potential");
        L->addWidget(new QLabel("Number of A_i coefficients:", D), 0, 0);
        L->addWidget(NAE, 0, 1);
        L->addWidget(new QLabel("Number of C_i coefficients:", D), 1, 0);
        L->addWidget(NPE, 1, 1);
        L->addWidget(ExcIntB, 2, 0, 1, 2);
        ExcIntB->setChecked(Worker->isExcIntAvailable());
        ExcIntB->setEnabled(Worker->isExcIntAvailable());
        L->setRowMinimumHeight(3, 20);
        L->addWidget(OK, 4, 0);
        L->addWidget(Cancel, 4, 1);
        connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
        connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
        do
        {
            if (D->exec() == QDialog::Rejected)
            {
                delete D;
                return 0.0;
            }
            NA = NAE->text().toInt();
            NC = NPE->text().toInt();
            if (NC < 3 || NC > 20) 
                QMessageBox::information(this, "MolSpektAnalysis", 
                    "For NA and NC you must choose valid integer numbers between 1 and 10 and 3 and 20, respectively");
        }
        while (NC < 3 || NC > 20 || NA <= 0 || NA > 10);
        ExcInt = ExcIntB->isChecked();
        delete D;
        getMinimum(Re, Te);
        getLRCoeffForReading(CNC, pCLRC, CLRC);
        if (Worker->isAnaPotAv() && ExcInt)
        {
            double A, b, g;
            AnaPot *aPot = Worker->getAnaPotForReading(); 
            aPot->getExchangeCoeff(A, b, g);
            Worker->unlock();
            TTPot->setExcInt(A, b, g);
        }
        TTPot->setNA(NA);
        if (CNC >= 3 ? !TTPot->calcPotCoeff(CLRC, 3, NC, Re, Te) : true)
        {
            delete TTPot;
            delete[] pCLRC;
            delete[] CLRC;
            QMessageBox::information(this, "MolSpektAnalysis", "Error calculating potential parameters!");
            return 0.0;
        }
        delete[] pCLRC;
        delete[] CLRC;
        double *R = new double[numPoints], *U = new double[numPoints], Uinf = Worker->getUinf(), Ri = getMinR(Uinf, 1e-3), Ra = 50.0;
        double *Sig = new double[numPoints];
        double h = (Ra - Ri) / (numPoints - 1), FQS;
        for (n=0; n < numPoints; n++)
        {
            R[n] = (n==0 ? Ri : R[n-1] + h);
            U[n] = getPoint(R[n]);
            Sig[n] = 1e-4;
        }
        nPot->setWorker(TTPot);
        State->addPotential(nPot);
        double b1 = 0.0, b2 = TTPot->getb(), b3 = b2 + 0.1;
        double FQS2 = TTPot->improvePotential(numPoints, R, U, Sig);
        double FQS1 = 0.0, FQS3 = TTPot->improvePotential(numPoints, R, U, Sig, b3);
        while (FQS3 < FQS2 && b3 < 4.5)
        {
            b1 = b2;
            FQS1 = FQS2;
            b2 = b3;
            FQS2 = FQS3;
            FQS3 = TTPot->improvePotential(numPoints, R, U, Sig, b3 += 0.1);
        }
        if (b1 == 0.0) 
        {
            FQS1 = TTPot->improvePotential(numPoints, R, U, Sig, b1 = b2 - 0.1);
            while (FQS1 < FQS2 && b1 >= 0.4)
            {
                b3 = b2;
                FQS3 = FQS2;
                b2 = b1;
                FQS2 = FQS1;
                FQS1 = TTPot->improvePotential(numPoints, R, U, Sig, b1 -= 0.1);
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
            FQS2 = TTPot->improvePotential(numPoints, R, U, Sig, b2 = b);
            ParabInterpol(b1, FQS1, b2, FQS2, b3, FQS3, b, FQS);
        }
        FQS = TTPot->improvePotential(numPoints, R, U, Sig, b);
        //b1 = 0.0; b2 = b; b3 = b2 + 0.01; 
        FQS2 = nPot->FitTangToenniesPot(false);
        /*TTPot->improvePotential(numPoints, R, U, Sig, b3);
        FQS3 = nPot->FitTangToenniesPot(false);
        while (FQS3 < FQS2 && b3 < 4.5)
        {
            b1 = b2;
            FQS1 = FQS2;
            b2 = b3;
            FQS2 = FQS3;
            TTPot->improvePotential(numPoints, R, U, Sig, b3 += 0.01);
            FQS3 = nPot->FitTangToenniesPot(false);
        }
        if (b1 == 0.0) 
        {
            TTPot->improvePotential(numPoints, R, U, Sig, b1 = b2 - 0.01);
            FQS1 = nPot->FitTangToenniesPot(false);
            while (FQS1 < FQS2 && b1 >= 0.4)
            {
                b3 = b2;
                FQS3 = FQS2;
                b2 = b1;
                FQS2 = FQS1;
                TTPot->improvePotential(numPoints, R, U, Sig, b1 -= 0.01);
                FQS1 = nPot->FitTangToenniesPot(false);
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
            TTPot->improvePotential(numPoints, R, U, Sig, b2 = b);
            FQS2 = nPot->FitTangToenniesPot(false);
            ParabInterpol(b1, FQS1, b2, FQS2, b3, FQS3, b, FQS);
        }
        TTPot->improvePotential(numPoints, R, U, Sig, b);
        FQS = nPot->FitTangToenniesPot();*/
        delete[] R;
        delete[] U;
        delete[] Sig;
        if (Worker->isAdCorrA())
        {
            int NAdCorr, TAdCorr, PAdCorr;
            double *adCorr, Rm, b, RIso1, RIso2;
            bool *adCorrFree;
            Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, &adCorrFree);
            nPot->Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, adCorrFree);
        }
        nPot->UpdateTab();
        nPot->setName(getName() + "TTPot");
        nPot->setSource("Potential from " + getSource());
        F = getFileName();
        nPot->setFileName(F.right(F.length() - F.indexOf('.') - 1) + "TTPot.pot");
        nPot->show();
        //return nPot->FitAnaPot();
        return FQS;
    }
    else
    {
        double FQS, WeightSum, *WF;
        int NWF, *WFR, NumFreePar;
        mWasMoving = false;
        if (molecule->getNumIso() == 0)
        {
            QMessageBox::information(this, "MolSpektAnalysis", "The potential has to be assigned to a molecule with atomic data first!");
            return 0.0;
        }
        setFitData();
        FQS = Worker->FitTangToenniesPot();
        if (showResult)
        {
            Worker->getFitResult(FQS, WeightSum, NWF, WFR, WF, NumFreePar);
            QMessageBox::information(this, tr("Fit result"), 
                                     "Initial chiSq=" + QString::number(FQSv, 'g', 6) + ", final chiSq="
                                     + QString::number(FQS, 'g', 6) + ", Sigma=" 
                                     + QString::number(sqrt(FQS/double(WeightSum - NumFreePar))));
            Worker->unlockBadList();
        }
        UpdateTab();
        Changed();
        return FQS;
    }
}

double Potential::FitIsoMass(int Iso, double InitStepSize, double& RIsoMass, int NumWFPoints)
{
    setFitData();
    return Worker->FitIsoMass(Iso, InitStepSize, RIsoMass, NumWFPoints);
}

void Potential::calcDerivatives(TableLine*, int, double*, double**)
{
    /*double h = (rmax - rmin) / double(NumWFPoints - 1), r, *Ev, **WF = 0, *F, *NF, S;
    int n, m, c, N, NC, NLRC, MaxPN = 0, J = -1, I = -1, Nv, *ri, *ra, *rc, v, NS = 0;
    if (LPot != 0) LPot->getNC(N, NC, NLRC);
    else aPot->getNC(N, NC, NLRC);
    double **DR = Create(NumWFPoints, N);
    if (LPot != 0) for (n=0, r = rmin; n < NumWFPoints; n++, r+=h)
    {
        LPot->getDerivatives(r, S, DR[n] + 3, DR[n] + NC + 3, DR[n][2], DR[n][1]);
        DR[n][0] = 1.0 - DR[n][1];
        DR[n][1] = -DR[n][0];
    }
    else for (n=0, r = rmin; n < NumWFPoints; n++, r+=h)
        aPot->getDerivatives(r, DR[n], DR[n] + NC);
    for (n=0; n < NL; n++)
    {
        if (Lines[n].Iso != I || Lines[n].Jss != J)
        {
            if (WF != 0)
            {
                delete[] Ev;
                Destroy(WF, Nv);
                delete[] F;
                delete[] NF;
                delete[] ri;
                delete[] ra;
                delete[] rc;
            }
            I = Lines[n].Iso;
            J = Lines[n].Jss;
            for (m=n, Nv = 1; (m < NL ? Lines[m].Iso == I && Lines[m].Jss == J : false); m++) if (Lines[m].vss >= Nv)
                Nv = Lines[m].vss + 1;
            getWaveFunc(I, J, Nv, Ev, WF, ri, ra, rc, F, NF);
        }
        if ((v = Lines[n].vss) >= Nv)
        {
            E[n] = -1.0;
            continue;
        }
        for (c=0; c < N; c++)
        {
            for (m = ri[v], S = 0.0; m <= rc[v]; m++) S += WF[v][m] * DR[m][c] * WF[v][m];
            for (S *= F[v] * F[v]; m < ra[v]; m++) S += WF[v][m] * DR[m][c] * WF[v][m];
            dyda[n][c] = NF[v] * S * NF[v];
        }
        E[n] = Ev[v];
        if (!Lines[n].isTE)
        {
            if (Lines[n].PN > MaxPN) MaxPN = Lines[n].PN;
            if (Lines[n].SourceN >= NS) NS = Lines[n].SourceN + 1;
        }
    }
    if (WF != 0)
    {
        delete[] Ev;
        Destroy(WF, Nv);
        delete[] F;
        delete[] NF;
        delete[] ri;
        delete[] ra;
        delete[] rc;
    }
    Destroy(DR, NumWFPoints);
    double W, **SE = Create(NS, MaxPN + 1), **WS = Create(NS, MaxPN + 1);
    //double ***SD = Create(NS, MaxPN + 1, N);
    for (m=0; m < NS; m++) for (n=0; n <= MaxPN; n++)
    {
        SE[m][n] = WS[m][n] = 0.0;
        //for (c=0; c<N; c++) SD[m][n][c] = 0.0;
    }
    for (n=0; n < NL; n++) if (!Lines[n].isTE)
    {
        E[n] *= -1.0;
        W = 1.0 / (Lines[n].err * Lines[n].err);
        SE[Lines[n].SourceN][Lines[n].PN] += W * E[n];
        WS[Lines[n].SourceN][Lines[n].PN] += W;
        //for (c=0; c < N; c++) SD[Lines[n].SourceN][Lines[n].PN][c] += W * dyda[n][c];
    }
    for (m=0; m < NS; m++) for (n=0; n <= MaxPN; n++) if (WS[m][n] != 0.0)
    {
        WS[m][n] = 1.0 / WS[m][n];
        SE[m][n] *= WS[m][n];
        //for (c=0; c<N; c++) SD[m][n][c] *= WS[m][n];
    }
    for (n=0; n < NL; n++) if (!Lines[n].isTE)
    {
        E[n] -= SE[Lines[n].SourceN][Lines[n].PN];
        for (c=0; c<N; c++) dyda[n][c] *= -1.0;// -= SD[Lines[n].SourceN][Lines[n].PN][c];
    }
    Destroy(SE, NS);
    Destroy(WS, NS);
    //Destroy(SD, NS, MaxPN + 1);*/
}

void Potential::setMolecule(Molecule *Mol)
{
    if (Fit->isRunning())
    {
        QMessageBox::information(this, "MolSpektAnalysis", "The molecular assignment cannot be changed while a fit is running!");
        return;
    }
    if (Mol != 0)
    {
        double UFact = 2e-18 * C_c * C_h * C_u / (C_hq * C_hq), *MU, *IsoF;
        int i, NIso;
        IsoTab *IsoT = Mol->getIso();
        if (IsoT != 0)
        {
            NIso = IsoT->numIso;
            MU = new double[NIso];
            IsoF = new double[NIso];
            for (i=0; i < NIso; i++)
            {    
                MU[i] = IsoT->redMass[i];
                IsoF[i] = MU[i] * UFact;
            }
            Worker->setMolData(IsoT, MU, IsoF);
        }
    }
    TableWindow::setMolecule(Mol);
}

double Potential::calcScatWaveFunc(bool write, double **Nodes, double maxR)        
{
    IsoTab *IsoT = Worker->getIsoT();
    if (IsoT == 0)
    {
        printf("Potential::calcScatWaveFunc(): Error: not enough data available!\n");
        return 0.0;
    }
    QTextStream S;
    if (write)
    {
        QString fileName = QFileDialog::getSaveFileName(this, tr("Write scattering wave functions"), 
                    MW->getDir(PotData), tr("Text tables (*.dat)"));
        QFile DDatei(fileName);
        if (!DDatei.open(QIODevice::WriteOnly))
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), tr("Error while opening file!"));
            return 0.0;
        }
        S.setDevice(&DDatei);
        S << "Scattering wave functions:\nR[A]";
    }
    int i, I, NI = IsoT->numIso, xi, Np;
    double RMin = 3.0, RMax = 103.0, R, M[NI], h = 0.001;
    double Y[NI][3], WF[NI][3], *U = getPoints(RMin, RMax, Np = int((RMax - RMin) / h) + 1);
    //printf("M=%f, h=%f, I=%f\n", M, h, IsoF[Iso]);
    bool gn[NI];
    int v[NI];
    double LE, AE, NE, *IsoF = Worker->getIsoF_ForReading();
    for (xi=0; U[xi] > 0.0 && xi < Np - 1; xi++) ;
    //printf("xi=%d\n", xi);
    for (LE = 0.0, AE = 1e-5, xi--; xi > 0 && AE < 1e10; xi--)
    {
        NE = AE * (2.0 + M[0] * U[xi]) - LE;
        LE = AE;
        AE = NE;
    }
    for (I=0; I < NI; I++) 
    {
        gn[I] = true;
        v[I] = 0;
        M[I] = IsoF[I] * h * h;
        Y[I][1] = 1e-5;
        Y[I][0] = Y[I][1] * exp(-sqrt(M[I]*U[xi]));
        WF[I][0] = Y[I][0] * PFakt(M[I], U[xi]);
        if (write) S << "\t" << IsoT->mNumIso1[I] << IsoT->chSymb1 << IsoT->mNumIso2[I] << IsoT->chSymb2;
    }
    //printf("xi=%d, U[xi]=%f\n", xi, U[xi]);
    if (write) S << "\n";
    for (i = xi + 1, R = RMin + h * double(xi); 
            R <= maxR && (RMin > 10.0 || i > xi + 50 ? fabs(M[0]*U[i]) > 1e-16 : true); i++)
    {
        //printf("Beginn\n");
        R += h;
        if (write) S << R;
        for (I=0; I < NI; I++)
        {
            //printf("i=%d\n", i);
            NCNF(M[I], U[i-1], Y[I] + 2, WF[I] + 2, -1);
            //printf("Nach NCNF\n");
            if (gn[I] && Y[I][2] < 0.0)
            {
                gn[I] = false;
                if (Nodes != 0 && v[I] <= cMaxv) Nodes[I][v[I]] = R - h * Y[I][2] / (Y[I][2] - Y[I][1]);
                v[I]++;
            }
            else if (!gn[I] && Y[I][2] > 0.0)
            {
                gn[I] = true;
                if (Nodes != 0 && v[I] <= cMaxv) Nodes[I][v[I]] = R - h * Y[I][2] / (Y[I][2] - Y[I][1]);
                v[I]++;
            }
            //printf("Vor write, WF[I]=%x, WF[I][2]=%f\n", WF[I], WF[I][2]);
            if (write) S << "\t" << QString::number(WF[I][1], 'e', 10);
            //printf("Nach write\n");
            Y[I][0] = Y[I][1];
            Y[I][1] = Y[I][2];
            WF[I][0] = WF[I][1];
            WF[I][1] = WF[I][2];
        }
        if (write) S << "\n";
        //if (i==18) break;
        if (i == Np - 1)
        {
            //printf("Vor delete\n");
            delete[] U;
            U = getPoints(RMin += 100.0, RMax += 100.0, Np);
            i=0;
        }
        //printf("Ende\n");
    }
    if (Nodes != 0) for (I=0; I < NI; I++) for (i = v[I]; i <= cMaxv; i++) Nodes[I][i] = 0.0;
    //for (I=0; I<NI; I++) 
        //printf("%d%s%d%s: Maxv=%d\n", IsoT->mNumIso1[I], IsoT->chSymb1.ascii(),
            //   IsoT->mNumIso2[I], IsoT->chSymb2.ascii(), v[I]);
    //printf("Vor delete\n");
    delete[] U;
    return R;
}

void Potential::NCNF(double M, double E, double *Y, double *W, int d)
{
    //printf("W=%x, Y=%x\n", W, Y);
    W[d] = Y[d] * PFakt(M, E);
    Y[0] = M * E * W[d] + 2.0 * Y[d] - Y[2*d];
}

double Potential::PFakt(double M, double E)
{
    return 1.0 / (1.0 - 0.0833333333333 * M * E);
}

bool Potential::readData(QString Filename)
{
    //printf("Potential::readData\n");
    if (CoupledComp != 0)
    {
        delete[] CoupledComp;
        CoupledComp = 0;
    }
    if (CoupledPot != 0)
    {
        delete[] CoupledPot;
        CoupledPot = 0;
    }
    Tab->setColumnCount(4);
    QFile Datei(Filename);
    if (!read(&Datei)) return false;
    Tab->blockSignals(true);
    int n, r, cc, lc = 0;
    bool Success = true, notSaved = false;
    QString Buffer;
    QTextStream S(&Datei);
    QString Spacer = " | ";
    QStringList L;
    mWasMoving = false;
    PreliminaryPotentialType potType = likelyASplinePotential;
    if ((Buffer = S.readLine()).left(8) != "Source: ") Success = false; 
    else 
    {
        setSource(Buffer.right(Buffer.length() - 8));
        if ((Buffer = S.readLine()).left(6) != "Name: ") Success = false;
        else
        {
            setName(Buffer.right(Buffer.length() - 6));
            for (r=0, cc = Tab->columnCount(); !S.atEnd(); r++)
            {
                if (Tab->rowCount() == r) Tab->setRowCount(r + 100);
                if ((Buffer = S.readLine()).left(15) == "Column titles: ")
                {
                    if (Buffer == "Column titles: R [A] | E [cm^-1] | d^2E/dR^2") potType = splinePot;
                    Buffer = S.readLine();
                }
                if (Buffer.left(21) == "Coupling information:") break;
                L = Buffer.split(Spacer);
                if ((lc = L.count()) > cc) Tab->setColumnCount(cc = lc);
                for (n=0; n < lc; n++) Tab->setItem(r, n, new QTableWidgetItem(L[n]));
                while (n < cc) Tab->setItem(r, n++, new QTableWidgetItem(""));
            }
            if (!S.atEnd())
            {
                r-=2;
                Buffer = S.readLine();
                QString MolFile, FB;
                if (Buffer.indexOf("CouplingFunction: ") == -1)
                {
                    Filename = Buffer.right(Buffer.length() - 16);
                    if (molecule != 0) Filename = getAbsolutePath(Filename, MolFile = molecule->getFileName());
                    S.readLine();
                    S.readLine();
                    CoupledComp = new int[10];
                    CoupledPot = new Potential*[10];
                    for (NCoupled = 0; NCoupled < 10 && !S.atEnd(); NCoupled++)
                    {
                        n = (Buffer = S.readLine()).indexOf(" | ");
                        if (n>0)
                        {
                            CoupledPot[NCoupled] = MW->getPotential((molecule != 0 ? getAbsolutePath(FB = Buffer.left(n), MolFile)
                                                                               : Buffer.left(n)), molecule);
                            if (CoupledPot[NCoupled] != 0 && CoupledPot[NCoupled]->didFileNameChange()) notSaved = true;
                                CoupledComp[NCoupled] = Buffer.right(Buffer.length() - n - 3).toInt();
                            if (CoupledPot[NCoupled] != 0 ? CoupledPot[NCoupled]->getCoupledSineWaveFunc() != 0 : false)
                                WaveFuncs = CoupledPot[NCoupled]->getCoupledSineWaveFunc();
                        }
                        else if (Buffer.left(18) == "CouplingFunction: ")
                        {
                            for (n=0; n < NCoupled; ++n) if (CoupledPot[n] != 0 && CoupledPot[n]->m_couplingFuncs != 0)
                                m_couplingFuncs = CoupledPot[n]->m_couplingFuncs;
                            if (m_couplingFuncs == 0)
                            {
                                Filename = Buffer.right(Buffer.length() - 18);
                                if (molecule != 0) Filename = getAbsolutePath(Filename, MolFile = molecule->getFileName());
                                m_couplingFuncs = new CouplingFuncs(MW);
                                m_couplingFuncs->readData(Filename);
                                MW->showMDIChild(m_couplingFuncs);
                            }
                            break;
                        }
                    }
                    if (WaveFuncs == 0)
                    {
                        WaveFuncs = new CoupledSineWaveFunc(MW, molecule);
                        if (!WaveFuncs->readData(Filename))
                        {
                            delete WaveFuncs;
                            WaveFuncs = 0;
                        }
                        else if (WaveFuncs->didFileNameChange()) notSaved = true;
                    }
                    if (WaveFuncs != 0) WaveFuncs->Assign();
                    if (!S.atEnd()) Buffer = S.readLine();
                }
                if (Buffer.left(18) == "CouplingFunction: ")
                {
                    Filename = Buffer.right(Buffer.length() - 18);
                    if (molecule != 0) Filename = getAbsolutePath(Filename, MolFile = molecule->getFileName());
                    if (CoupledPot == 0)
                    {
                        CoupledPot = new Potential*[10];
                        while (!S.atEnd())
                        {
                            Buffer = S.readLine();
                            if (Buffer.left(19) == "Coupled potential: ")
                            {
                                FB = Buffer.right(Buffer.length() - 19);
                                CoupledPot[NCoupled++] = MW->getPotential((molecule != 0 ? getAbsolutePath(FB, MolFile) : FB), molecule);
                            }
                        }
                    }
                    for (n=0; m_couplingFuncs != 0 && n < NCoupled; ++n) if (CoupledPot[n] != 0 && CoupledPot[n]->m_couplingFuncs != 0)
                        m_couplingFuncs = CoupledPot[n]->m_couplingFuncs;
                    if (m_couplingFuncs == 0)
                    {
                        m_couplingFuncs = new CouplingFuncs(MW);
                        m_couplingFuncs->readData(Filename);
                        MW->showMDIChild(m_couplingFuncs);
                    }
                }
            }
            if (lc == 0) r--;
            Tab->setRowCount(r);
        }
    }
    if (!Success) 
    {
        Tab->blockSignals(false);
        if (!readPointPot(Filename)) if (!readPotData(Filename)) if (!readCoupledPotData(Filename)) return false;
    }
    else 
    {
        //printf("TableWindow::readData successfull\n");
        UpdatePot(potType);
    }
    Tab->blockSignals(false);
    if (notSaved) Changed();
    else Saved();
    //printf("Ende Potential::readData\n");
    return true;
}

bool Potential::testIfSplinePot()
{
    bool E1, E2;
    double D1 = Tab->item(0, 0)->text().toDouble(&E1), D2 = Tab->item(1, 0)->text().toDouble(&E2);
    return E1 && E2 && D1 < D2;
}

double **Potential::fitNodePos()
{
    AnaPot *aPot = Worker->getAnaPotForWriting();
    if (aPot == 0) 
    {
        QMessageBox::information(this, tr("QT4MolSpektAn"), 
                    tr("fitNodePos: fatal error, no analytical potential data available or a different fit is running!"));
        return 0;
    }
    IsoTab *IsoT = Worker->getIsoT();
    if (IsoT == 0)
    {
        QMessageBox::information(this, tr("QT4MolSpektAn"),
                                 tr("fitNodePos: fatal error, no isotopic data available!"));
        return 0;
    }
    int NI = IsoT->numIso, NC = Tab->rowCount();
    double D, P, oN, R, F;
    double b, De, iCoeff, iOff, iExp, Ra, Ri, Rm, Tm, *LRCoeff, *PotCoeff;
    double **NodePos = Create(NI, cMaxv+1);
    int nPotCoeff, nLRCoeff, *pLRCoeff, n;
    bool error = true;
    aPot->getPotCoeff(b, Rm, Tm, De, nPotCoeff, PotCoeff);
    aPot->getInnerWall(Ri, iOff, iExp, iCoeff);
    aPot->getLRCoeffForWriting(nLRCoeff, pLRCoeff, LRCoeff, 0, &Ra);
    QDialog *NFDialog = new QDialog(this);
    NFDialog->setMinimumSize(300, 270);
    NFDialog->setMaximumSize(300, 270);
    QLabel *Label1 = new QLabel("Isotopomer:", NFDialog);
    Label1->setGeometry(10, 10, 120, 20);
    QLabel *Label2 = new QLabel("Number of node:", NFDialog);
    Label2->setGeometry(10, 40, 120, 20);
    QLabel *Label3 = new QLabel("Value to fit to [A]:", NFDialog);
    Label3->setGeometry(10, 70, 120, 20);
    QLabel *Label6 = new QLabel("Parameter to vary:", NFDialog);
    Label6->setGeometry(10, 100, 120, 20);
    QLabel *Label4 = new QLabel("Wanted accuracy [A]:", NFDialog);
    Label4->setGeometry(10, 130, 120, 20);
    QLabel *Label5 = new QLabel("Max search radius [A]:", NFDialog);
    Label5->setGeometry(10, 160, 120, 20);
    QLabel *Label7 = new QLabel("Max iterations:", NFDialog);
    Label7->setGeometry(10, 190, 120, 20);
    QComboBox *isotopomer = new QComboBox(NFDialog);
    isotopomer->setGeometry(140, 10, 150, 20);
    QLineEdit *number = new QLineEdit(QString::number(fitN_N), NFDialog);
    number->setGeometry(140, 40, 150, 20);
    QLineEdit *value = new QLineEdit(QString::number(fitN_V, 'f', 1), NFDialog);
    value->setGeometry(140, 70, 150, 20);
    QComboBox *parameter = new QComboBox(NFDialog);
    parameter->setGeometry(140, 100, 150, 20);
    QLineEdit *accuracy = new QLineEdit(QString::number(fitN_A, 'f', 3), NFDialog);
    accuracy->setGeometry(140, 130, 150, 20);
    QLineEdit *maxr = new QLineEdit(QString::number(fitN_maxr, 'f', 1), NFDialog);
    maxr->setGeometry(140, 160, 150, 20);
    QLineEdit *maxIt = new QLineEdit(QString::number(fitN_MI), NFDialog);
    maxIt->setGeometry(140, 190, 150, 20);
    QPushButton *OK = new QPushButton("OK", NFDialog);
    OK->setGeometry(10, 240, 80, 20);
    QPushButton *Cancel = new QPushButton("Cancel", NFDialog);
    Cancel->setGeometry(210, 240, 80, 20);
    for (n=0; n < NI; n++) 
        isotopomer->addItem(QString::number(IsoT->mNumIso1[n]) + *IsoT->chSymb1 
                + QString::number(IsoT->mNumIso2[n]) + *IsoT->chSymb2);
    isotopomer->setEditable(false);
    isotopomer->setCurrentIndex(fitN_I);
    for (n=0; n < nPotCoeff; n++) parameter->addItem(Tab->item(n, 0)->text());
    parameter->addItem(Tab->item(nPotCoeff + 3, 0)->text());
    for (n = nPotCoeff + 6; n < NC - 3; n++) parameter->addItem(Tab->item(n, 0)->text());
    parameter->setEditable(false);
    parameter->setCurrentIndex(fitN_C);
    connect(OK, SIGNAL(clicked()), NFDialog, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), NFDialog, SLOT(reject()));
    while (error)
    {
        if (NFDialog->exec() == QDialog::Rejected) return 0;
        fitN_I = isotopomer->currentIndex();
        fitN_N = number->text().toInt();
        if (fitN_N <= 0 || fitN_N > cMaxv) 
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                     tr("The selected number of the node is invalid!"));
            continue;
        }
        if ((fitN_V = value->text().toDouble()) <= rmin)
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                     tr("The selected value to fit to is invalid!"));
            continue;
        }
        fitN_C = parameter->currentIndex();
        if ((fitN_A = accuracy->text().toDouble()) < 1e-15)
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                     tr("The selected accuracy is invalid!"));
            continue;
        }
        if ((fitN_maxr = maxr->text().toDouble()) < fitN_V)
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                tr("The selected maximum radius is smaller than the value to fit to!"));
            continue;
        }
        fitN_MI = maxIt->text().toInt();
        if ((R = calcScatWaveFunc(false, NodePos, fitN_maxr)) < fitN_V)
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                    QString("The selected value to fit to is outside the R range which can be calculated of Rmax=%1!").arg(R, 6, 'f', 4));
            continue;
        }
        if (NodePos[fitN_I][fitN_N-1] == 0.0) 
        {
            QMessageBox::information(this, tr("QT4MolSpektAn"), tr("The selected node does not exist!"));
            continue;
        }
        error = false;
    }
    if (fitN_C < nPotCoeff) P = PotCoeff[fitN_C];
    else if (fitN_C == nPotCoeff) P = Ri;
    else if (fitN_C == nPotCoeff + 1) P = iCoeff;
    else if (fitN_C == nPotCoeff + 2) P = Ra;
    else P = LRCoeff[fitN_C - nPotCoeff - 3];
    P += (D = 0.1 * P);
    printf("fitN_MI=%d, Node=%f, fitN_V=%f, fitN_A=%f\n", fitN_MI, NodePos[fitN_I][fitN_N-1], fitN_V, fitN_A);
    for (n=0; n < fitN_MI && fabs((oN = NodePos[fitN_I][fitN_N-1]) - fitN_V) > fitN_A; n++)
    {
        printf("n=%d, D=%g, P=%g\n", n, D, P);
        printf("Node=%g, V=%g\n", NodePos[fitN_I][fitN_N-1], fitN_V);
        if (fitN_C < nPotCoeff) PotCoeff[fitN_C] = P;
        else if (fitN_C == nPotCoeff) 
        {
            Ri = P;
            aPot->setInnerWall(Ri, iOff, iExp, iCoeff);
        }
        else if (fitN_C == nPotCoeff + 1) 
        {
            iCoeff = P;
            aPot->setInnerWall(Ri, iOff, iExp, iCoeff);
        }        
        else if (fitN_C == nPotCoeff + 2) 
        {
            Ra = P;
            aPot->setLRCoeff(Ra, nLRCoeff, pLRCoeff, LRCoeff);
        }
        else LRCoeff[fitN_C - nPotCoeff - 3] = P;
        aPot->sConnect();
        calcScatWaveFunc(false, NodePos, fitN_maxr);
        if (NodePos[fitN_I][fitN_N-1] == 0.0) D -= (P *= 0.5);
        else 
        {
            if ((F = (fitN_V - NodePos[fitN_I][fitN_N-1]) / (NodePos[fitN_I][fitN_N-1] - oN)) > 10.0) F = 10.0;
            else if (F < -10.0) F = -10.0;
            P += (D *= F);
        }
    }
    if (fabs((oN = NodePos[fitN_I][fitN_N-1]) - fitN_V) > fitN_A) 
        QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                 tr("The fit was not successful within the given number of iterations!"));
    else QMessageBox::information(this, tr("QT4MolSpektAn"), 
                                 tr("The fit was successful!"));
    UpdateTab();
    delete NFDialog;
    Changed();
    return NodePos;
}

void Potential::getTexTable(int NumWFPoints, double FQS)
{
    double b, De, iCoeff, iOff, iExp, Ra, Ri, Rm, Tm, *LRCoeff, *PotCoeff, MinR, MinU, U_inf;
    int nPotCoeff, nLRCoeff, *pLRCoeff, n, m, l, p, q, prec = 0, lang, nDig, r;
    QStringList L;
    QString Buffer, eB, B2;
    QTableWidgetItem *I;
    QDialog *D = new QDialog(this);
    QGridLayout *La = new QGridLayout(D);
    D->setWindowTitle("Select language and precision");
    La->addWidget(new QLabel("Precision [cm^-1]:", D), 0, 0);
    QLineEdit *PrecE = new QLineEdit("0.0001", D);
    La->addWidget(PrecE, 0, 1);
    La->addWidget(new QLabel("Language:", D), 1, 0);
    QComboBox *LangB = new QComboBox(D);
    LangB->setEditable(false);
    LangB->addItems(QStringList() << "english" << "german");
    La->addWidget(LangB, 1, 1);
    La->setRowMinimumHeight(2, 20);
    QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
    La->addWidget(OK, 3, 0);
    La->addWidget(Cancel, 3, 1);
    connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
    if (D->exec() == QDialog::Rejected)
    {
        delete D;
        return;
    }
    double Prec = PrecE->text().toDouble(), x, X;
    if (Prec > 0.0)
    {
        Prec = 1.0 / Prec;
        nDig = int(log10(Prec));
        if (pow(10, nDig) != Prec) nDig++;
    }
    else nDig = -1;
    lang = LangB->currentIndex();
    delete D;
    L << "\\begin{table}" << "\\caption{???}";
    L << "\\label{" + (State != 0 ? State->getName() : QString("???"))  + "Potential}";
    L << "\\centering";
    Tab->blockSignals(true);
    if (Worker->isAnaPotAv())
    {
        AnaPot *aPot = Worker->getAnaPotForReading();
        double q_e, q_f;
        int intIExp;
        Worker->getLambdaDoubling(q_e, q_f);
        aPot->getPotCoeff(b, Rm, Tm, De, nPotCoeff, PotCoeff);
        aPot->getInnerWall(Ri, iOff, iExp, iCoeff);
        aPot->getLRCoeffForReading(nLRCoeff, pLRCoeff, LRCoeff, &Ra, true);
        if ((x = abs((Ri - Rm) / (Ri + b * Rm))) < (X = abs((Ra - Rm) / (Ra + b * Rm)))) x=X;
        L << "\\begin{tabular}{lp{0.1\\textwidth}r}\\hline\\noalign{\\smallskip}";
        for (n=0, X=x; n < nPotCoeff; n++, X*=x) 
        {
            m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
            if (Prec > 0.0 && m<0)
            {
                prec = int(log10(Prec * X * abs(PotCoeff[n])));
                if (prec < 0) prec = 0;
            }
            m = (Buffer = QString::number(PotCoeff[n], 'e', (m>=0 ? m : (Prec > 0.0 ? prec : 16)))).indexOf("e");
            Tab->item(n, 1)->setText(Buffer);
            Buffer = "$a_{" + QString::number(n+1) + "}$   && " + Buffer.left(m) + "$\\times 10^{" 
                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) + "}$ cm$^{-1}$ \\\\";
            L << (lang == 1 ? Buffer.replace('.', ',') : Buffer);
        }
        L << "\\hline\\noalign{\\smallskip}";
        Buffer = QString::number(b, 'g');
        if (Buffer.indexOf('.') == -1) Buffer += ".0";
        Tab->item(n++, 1)->setText(Buffer);
        L << "$b$   && " 
                + (lang == 1 ? Buffer.replace('.', ',') : Buffer)
                + "    \\\\";
        Buffer = QString::number(Rm, 'g', 12);
        if (Buffer.indexOf('.') == -1) Buffer += ".0";
        Tab->item(n++, 1)->setText(Buffer);
        L << "$R_m$ && " 
                + (lang == 1 ? Buffer.replace('.', ',') : Buffer)
                + " \\AA \\\\ ";
        Buffer = QString::number(Tm, 'f', ((m = (I = Tab->item(n, 3)) != 0 ? I->text().toInt() : 0) > 0 ? m 
                    : (nDig >= 0 ? nDig : 4)));
        Tab->item(n++, 1)->setText(Buffer);
        L << "$T_m$ && " 
                + (lang == 1 ? Buffer.replace('.', ',') : Buffer)
                + " cm$^{-1}$ \\\\ \\hline\\noalign{\\smallskip}";
        Buffer = QString::number(Ri, 'g');
        if (Buffer.indexOf('.') == -1) Buffer += ".0";
        Tab->item(n++, 1)->setText(Buffer);
        L << "$R_i$ && " 
                + (lang == 1 ? Buffer.replace('.', ',') : Buffer)
                + " \\AA \\\\";
        m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() : 0);
        if (Prec > 0.0 && m<=0)
        {
            prec = int(log10(Prec * iExp * iExp * iCoeff * pow(Ri, -iExp - 1.0)));
            if (prec < 0) prec = 0;
        }
        Buffer = QString::number(iExp, 'g', (m>0 ? m : (Prec > 0.0 ? prec + 1 : 11)));
        Tab->item(n++, 1)->setText(Buffer);
        intIExp = Buffer.toInt();
        L << "$s$   && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + "    \\\\";
        m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() : -1);
        if (Prec > 0.0 && m<=0)
        {
            prec = int(log10(Prec * abs(iOff)));
            if (prec < 0) prec = 0;
        }
        m = (Buffer = QString::number(iOff, 'e', (m>0 ? m : (Prec > 0.0 ? prec : 4)))).indexOf("e");
        Tab->item(n++, 1)->setText(Buffer);
        L << "$A$   && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) 
                + "}$ cm$^{-1}$ \\\\";
        m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() : -1);
        if (Prec > 0.0 && m<=0)
        {
            prec = int(log10(abs(Prec * iCoeff * pow(Ri, -iExp))));
            if (prec < 0) prec = 0;
        }
        m = (Buffer = QString::number(iCoeff, 'e', (m>0 ? m : (Prec > 0.0 ? prec : 4)))).indexOf("e");
        Tab->item(n++, 1)->setText(Buffer);
        L << "$B$   && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) 
                + "}$ cm$^{-1}$\\AA$^" + (intIExp > 0 ? QString::number(intIExp) : "s") 
                + "$ \\\\ \\hline\\noalign{\\smallskip}";
        Buffer = QString::number(Ra, 'g');
        Tab->item(n++, 1)->setText(Buffer);
        L << "$R_a$ && " 
                + (lang == 1 ? Buffer.replace('.', ',') : Buffer) 
                + " \\AA \\\\";
        for (l=0; l < nLRCoeff; l++)
        {
            m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
            if (Prec > 0.0 && m<0)
            {
                prec = int(log10(abs(Prec * LRCoeff[l] * pow(Ra, -pLRCoeff[l]))));
                if (prec < 0) prec = 0;
            }
            m = (Buffer = QString::number(LRCoeff[l], 'e', (m>=0 ? m : (Prec > 0.0 ? prec : 16)))).indexOf("e");
            Tab->item(n++, 1)->setText(Buffer);
            eB = QString::number(pLRCoeff[l]);
            L << "$C_{" + eB + "}$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                    + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) + "}$ cm$^{-1}$\\AA$^{" 
                    + eB + "}$ \\\\";
        }
        //if (De != 0.0) 
        //{
        Buffer = QString::number(De, 'f', ((m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() : -1)) > 0 ? m 
                            : (Prec > 0 ? nDig : 4)));
        if (Buffer.indexOf('.') == -1) Buffer += ".0";
        Tab->item(n++, 1)->setText(Buffer);
        L << "$U_\\infty$ && " 
            + (lang == 1 ? Buffer.replace('.', ',') : Buffer)
            + " cm$^{-1}$ \\\\ \\hline\\noalign{\\smallskip}";
        //}
        /*L << "R$_{min}$   && " + QString::number(rmin, 'f', 
                                    ((I = Tab->item(++n, 3)) != 0 ? I->text().toInt() : 1)) + " \\AA \\\\";
        L << "R$_{max}$   && " 
                + QString::number(rmax, 'f', ((I = Tab->item(++n, 3)) != 0 ? I->text().toInt(): 1)) 
                + " \\AA \\\\ \\hline";*/
        if (q_e != 0.0 || q_f != 0.0)
        {
            double *IsoF = Worker->getIsoF_ForReading();
            l = (State != 0 ? (State->getFitData() != 0 ? State->getFitData()->getMaxJ() : 300) : 300);
            if (q_e != 0.0)
            {
                m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
                if (Prec > 0.0 && m<0)
                {
                    prec = int(log10(abs(Prec * q_e * double(l * (l+1)) / (Ri * Ri * IsoF[0])))) + 1;
                    if (prec < 0) prec = 0;
                }
                m = (Buffer = QString::number(q_e, 'g', (m>=0 ? m+1 : (Prec > 0.0 ? prec : 5)))).indexOf("e");
                if (Buffer.indexOf('.') == -1) Buffer += ".0";
                Tab->item(n++, 1)->setText(Buffer);
                if (m>0)  L << "$\\eta_e$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) + "}$\\\\";
                else L << "$\\eta_e$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + "\\\\";
            }
            if (q_f != 0.0)
            {
                m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
                if (Prec > 0.0 && m<0)
                {
                    prec = int(log10(abs(Prec * q_f * double(l * (l+1)) / (Ri * Ri * IsoF[0])))) + 1;
                    if (prec < 0) prec = 0;
                }
                m = (Buffer = QString::number(q_f, 'g', (m>=0 ? m+1 : (Prec > 0.0 ? prec : 5)))).indexOf("e");
                if (Buffer.indexOf('.') == -1) Buffer += ".0";
                Tab->item(n++, 1)->setText(Buffer);
                if (m>0) L << "$\\eta_f$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) + "}$\\\\";
                else L << "$\\eta_f$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + "\\\\";
            }
            L << "\\hline\\noalign{\\smallskip}";
        }
        else l=-1;
        r=n;
        if (Worker->isAdCorrA())
        {
            int NAdCorr, TAdCorr, PAdCorr;
            double *adCorr, Rm, b, RIso1, RIso2;
            Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, 0, true);
            X = pow(2.0 * Rm / (Ri + Rm), PAdCorr);
            if (NAdCorr > 1) x = (Ri - Rm) / (Ri + b * Rm);
            for (n=0; n < NAdCorr; n++, X*=x)
            {
                m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
                if (Prec > 0.0 && m<0)
                {
                    prec = int(log10(abs(Prec * adCorr[n] * X))) + 1;
                    if (prec < 0) prec = 0;
                }
                m = (Buffer = QString::number(adCorr[n], 'g', (m>=0 ? m+1 : (Prec > 0.0 ? prec : 10)))).indexOf('e');
                if (Buffer.indexOf('.') == -1) Buffer += ".0";
                Tab->item(r++, 1)->setText(Buffer);
                if (m>0) L << "$v_{" + QString::number(n) + "}$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) 
                                + "$\\times 10^{" + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) 
                                    + "}$ cm$^{-1}$\\\\";
                else L << "$v_{" + QString::number(n) + "}$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) 
                            + " cm$^{-1}$\\\\"; 
            }
            L << "\\hline\\noalign{\\smallskip}";
            r+=4;
            delete[] adCorr;
        }
        if (Worker->getNSpinRGamma() > 0)
        {
            int NSpinRGamma;
            double *SpinRGamma, SpinRR1, SpinRR2, SpinRb, SpinRRm, *IsoF = Worker->getIsoF_ForReading();
            Worker->getSpinRGamma(NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2, SpinRb, SpinRRm);
            if (l==-1) l = (State != 0 ? (State->getFitData() != 0 ? State->getFitData()->getMaxJ() : 300) : 300);
            X = Prec * double(l+1) / (Ri * Ri * IsoF[0]);
            if (NSpinRGamma > 1) x = (Ri - Rm) / (Ri + b * Rm);
            for (n=0; n < NSpinRGamma; n++, X*=x)
            {
                m = ((I = Tab->item(n, 3)) != 0 ? I->text().toInt() - 1 : -1);
                if (Prec > 0.0 && m<0)
                {
                    prec = int(log10(abs(X * SpinRGamma[n]))) + 1;
                    if (prec < 0) prec = 0;
                }
                m = (Buffer = QString::number(SpinRGamma[n], 'g', (m>=0 ? m+1 : (Prec > 0.0 ? prec : 6)))).indexOf('e');
                if (Buffer.indexOf('.') == -1) Buffer += ".0";
                Tab->item(r++, 1)->setText(Buffer);
                if (m>0) L << "$\\gamma_{" + QString::number(n) + "}$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m)
                                + "$\\times 10^{" + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) + "}$\\\\";
                else L << "$\\gamma_{" + QString::number(n) + "}$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + "\\\\";
            }
            Buffer = QString::number(SpinRR1, 'g');
            if (Buffer.indexOf('.') == -1) Buffer += ".0";
            L << "$R_1$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + " \\AA \\\\";
            Buffer = QString::number(SpinRR2, 'g');
            if (Buffer.indexOf('.') == -1) Buffer += ".0";
            L << "$R_2$ && " + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + " \\AA \\\\";
            L << "\\hline\\noalign{\\smallskip}";
        }
        L << (lang == 0 ? "\\multicolumn{2}{l}{derived constants:}       \\\\" 
                        : "\\multicolumn{2}{l}{Abgeleitete Konstanten:}  \\\\");
        Worker->unlock();
        aPot->getMinimum(MinR, MinU);
        Buffer = (De == 0.0 ? "$D_e$    && " + QString::number(-1.0 * MinU, 'f', 3) :
                "$T_e$    && " + QString::number(MinU, 'f', 4)) + " cm$^{-1}$ \\\\";
        L << (lang == 1 ? Buffer.replace('.', ',') : Buffer);
        L << "$R_e$    && " + (lang == 1 ? QString::number(MinR, 'f', 5).replace('.', ',') 
                                : QString::number(MinR, 'f', 5)) + " \\AA \\\\ \\hline";
    }
    else if (Worker->isLPot())
    {
        MLRPot *LPot = dynamic_cast<MLRPot*>(Worker);
        LPot->getCoefficients(p, q, U_inf, De, MinR, Rm, nPotCoeff, PotCoeff, nLRCoeff, 
                              pLRCoeff, LRCoeff);
        L << "\\begin{tabular}{lp{0.1\\textwidth}r}\\hline\\noalign{\\smallskip}";
        for (n=0; n < nPotCoeff; n++)
            L << ("$\\beta_{" + QString::number(n) + "}$   ").left(14) + "&&"
                 + ("  " + QString::number(PotCoeff[n], 'f', 12)).right(17) + " \\\\";
        for (n=0; n < nLRCoeff; n++)
            L << ("$C_{" + QString::number(pLRCoeff[n]) + "}$       ").left(14) + "&& "
                 + QString::number(LRCoeff[n], 'e', 3).replace("e", "$\\times 10^{") 
                 + "}$ cm$^{-1}$\\AA$^{" + QString::number(pLRCoeff[n]) + "}$ \\\\";
        if (U_inf > 0.0) 
            L << "$U_\\infty     && " + QString::number(U_inf, 'f', 4) + " cm$^{-1}$ \\\\"; 
        L << "$D_e$         && " + QString::number(De, 'f', 4) + " cm$^{-1}$ \\\\";
        L << "$R_e$         && " + QString::number(MinR, 'f', 7) + " \\AA \\\\";
        L << "$R_{ref}$     && " + QString::number(Rm, 'f', 1) + " \\AA \\\\";
        L << "$p$           && " + QString::number(p) + " \\\\"; 
        L << "$q$           && " + QString::number(q) + " \\\\ \\hline";
    }
    else if (Worker->isSplinePot())
    {
        int numSplinePoints;
        if (!getTexTableSplineFitRunning)
        {
            FQSv = getFQS(NumWFPoints);
            SplinePoint *points = Worker->getSplinePoints(numSplinePoints);
            for (n=0; n < numSplinePoints; n++) 
                Tab->item(n, 0)->setText(QString::number(points[n].x, 'f', 2));
            delete points;
            UpdatePot(splinePot);
            getTexTableSplineFitRunning = true;
            FitSplinePot(false);
            return;
        }
        if (QMessageBox::Yes == QMessageBox::question(this, "MolSpektAnalysis", "The new chisq is " + QString::number(FQS)
                + " while the old chisq was " + QString::number(FQSv) + ", do you want one more fit iteration?", 
                QMessageBox::Yes | QMessageBox::No))
        {
            FitSplinePot(false);
            return;
        }
        else getTexTableSplineFitRunning = false;
        int NLRC, *pLRC;
        SplinePoint *points;
        double *LRC, iA, iO, iExp, Uinf = Worker->getUinf();
        Worker->getSplinePotForReading(numSplinePoints, points, NLRC, pLRC, LRC, iA, iO, iExp);
        if ((m = numSplinePoints / 2) * 2 < numSplinePoints) m++;
        L << "\\begin{tabular}{ll|ll}";
        L << "\\hline\\noalign{\\smallskip}";
        if (lang == 0)
            L << "R [\\AA] & Energy [cm$^{-1}$] & R [\\AA] & Energy [cm$^{-1}$]\\\\";
        else L << "R [\\AA] & Energie [cm$^{-1}$] & R [\\AA] & Energie [cm$^{-1}$]\\\\";
        L << "\\hline\\noalign{\\smallskip}";
        if (lang == 0) for (n=0; n<m; n++)
        {
            Buffer = QString::number(points[n].y, 'f', nDig);
            Tab->item(n, 1)->setText(Buffer);
            L << QString::number(points[n].x, 'f', 2) + " & " + Buffer;
            Buffer = QString::number(points[n+m].y, 'f', nDig);
            Tab->item(n+m, 1)->setText(Buffer);
            L << (m+n < numSplinePoints ? " & " + QString::number(points[n+m].x, 'f', 2) + " & "
                                            + Buffer + " \\\\" 
                                        : " &        &           \\\\");
        }
        else for (n=0; n<m; n++)
        {
            Buffer = QString::number(points[n].y, 'f', nDig);
            Tab->item(n, 1)->setText(Buffer);
            L << QString::number(points[n].x, 'f', 2).replace('.', ',') + " & " 
                    + Buffer.replace('.', ',');
            Buffer = QString::number(points[n+m].y, 'f', nDig);
            Tab->item(n+m, 1)->setText(Buffer);
            L << (m+n < numSplinePoints ? 
                      " & " + QString::number(points[n+m].x, 'f', 2).replace('.', ',') 
                      + " & " + Buffer.replace('.', ',') + " \\\\" 
                    : " &        &           \\\\");
        }
        L << "\\hline";
        L << "\\end{tabular}";
        L << "\\begin{tabular}{ll|ll} \\hline\\noalign{\\smallskip}";
        l = NLRC + (Uinf > 0.0 ? 4 : 3);
        QString Buff[l];
        if (Prec > 0.0)
        {
            prec = int(log10(abs(Prec * iO)));
            if (prec < 0) prec = 0;
        }
        m = (Buffer = QString::number(iO, 'e', (Prec > 0.0 ? prec : 4))).indexOf("e");
        Tab->item(numSplinePoints + 1, 1)->setText(Buffer);
        Buff[0] = "$A$ & " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) 
                + "}$ cm$^{-1}$";
        if (Prec > 0.0)
        {
            prec = int(log10(Prec * iA * pow(points[0].x, -6.0)));
            if (prec < 0) prec = 0;
        }
        m = (Buffer = QString::number(iA, 'e', (Prec > 0.0 ? prec : 4))).indexOf("e");
        Tab->item(numSplinePoints, 1)->setText(Buffer);
        Buff[1] = "$B$ & " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(m) + "$\\times 10^{" 
                + QString::number(Buffer.right(Buffer.length() - m - 1).toInt()) 
                + "}$ cm$^{-1}$\\AA$^6$";
        Buff[n=2] = "$n$ & 6";
        if (Uinf > 0.0)
        {
            Buffer = QString::number(Uinf, 'f', (Prec > 0 ? nDig : 3));
            Tab->item(numSplinePoints + n, 1)->setText(Buffer);
            Buff[++n] = "$U_\\infty $ & " 
                         + (lang == 1 ? Buffer.replace('.', ',') : Buffer) + " cm$^{-1}$";
        }
        for (m=0, Ra = points[numSplinePoints - 1].x; m < NLRC; m++)
        {
            if (Prec > 0.0)
            {
                prec = int(log10(Prec * abs(LRC[m]) * pow(Ra, -pLRC[m])));
                if (prec < 0) prec = 0;
            }
            p = (Buffer = QString::number(LRC[m], 'e', (Prec > 0.0 ? prec : 16))).indexOf("e");
            Tab->item(numSplinePoints + n, 1)->setText(Buffer);
            eB = QString::number(pLRC[m]);
            Buff[++n] = "$C_{" + eB + "}$ & " + (lang == 1 ? Buffer.replace('.', ',') : Buffer).left(p) + "$\\times 10^{" 
                      + QString::number(Buffer.right(Buffer.length() - p - 1).toInt()) + "}$ cm$^{-1}$\\AA$^{" 
                      + eB + "}$";
        }
        if ((m = l / 2) * 2 < l) m++;
        for (n=0; n<m; n++) L << Buff[n] + (n + m < l ? " & " + Buff[n+m] + " \\\\"
                                                      : " &       &     \\\\");
        L << "\\hline";
        Worker->unlock();
    }
    L << "\\end{tabular}" << "\\end{table}";
    updatePot();
    Tab->blockSignals(false);
    Changed();
    QTextEdit *Window = new QTextEdit();
    Window->setPlainText(L.join("\n"));
    MW->showMDIChild(Window);
}

bool Potential::writeData(QString Filename, bool writeFitData)
{
    m_saving = true;
    QFile Datei(Filename);
    if (!write(&Datei)) return false;
    int r, c, R = Tab->rowCount(), C = Tab->columnCount();
    QTextStream S(&Datei);
    QStringList L;
    QTableWidgetItem *I;
    QString Spacer = " | ";
    S << "Source: " << getSource() << "\n";
    S << "Name: " << getName() << "\n";
    for (c=0; c<C; c++) L << ((I = Tab->horizontalHeaderItem(c)) != 0 ? I->text() : "?");
    S << "Column titles: " << L.join(Spacer) << "\n";
    for (r=0; r < R; r++) 
    {
        L.clear();
        for (c=0; c < C; c++) L << ((I = Tab->item(r, c)) != 0 ? I->text() : "");
        S << L.join(Spacer) << "\n";
    }
    if (WaveFuncs != 0)
    {
        if (!WaveFuncs->isSaved())
        {
            if (WaveFuncs->getFileName().isEmpty())
            {
                Filename = WaveFuncs->getSource();
                Filename = Filename.right(Filename.length() - Filename.lastIndexOf(QRegExp("[\\/]")) - 1);
                Filename = Filename.left(Filename.indexOf('.') + 1) + "SWF";
                WaveFuncs->setFileName(Filename);
            }
            WaveFuncs->writeData();
        }
        if (NCoupled > 0 && WaveFuncs->isSaved())
        {
            QString MolFile;
            S << "\nCoupling information:\n";
            S << "Wave functions: " << (molecule != 0 ? WaveFuncs->getRelativePath(MolFile = molecule->getFileName())
                                                      : WaveFuncs->getFileName()) << '\n';
            S << "\nCoupled potential: | Coupled component:\n";
            for (r=0; r < NCoupled; r++) if (CoupledPot[r] != 0)
            {
                if (!CoupledPot[r]->isSaved() && !CoupledPot[r]->Saving()) CoupledPot[r]->writeData();
                S << (molecule != 0 ? CoupledPot[r]->getRelativePath(MolFile) : CoupledPot[r]->getFileName(true))
                  << " | " << QString::number(CoupledComp[r]) << '\n';
            }
        }
    }
    if (m_couplingFuncs != 0)
    {
        if (!m_couplingFuncs->isSaved()) m_couplingFuncs->writeData();
        QString MolFile;
        S << "CouplingFunction: " << (molecule != 0 ? m_couplingFuncs->getRelativePath( MolFile = molecule->getFileName()) : m_couplingFuncs->getFileName()) << '\n';
        if (WaveFuncs == 0) for (r=0; r < NCoupled; ++r) if (CoupledPot[r] != 0)
        {
            if (!CoupledPot[r]->isSaved() && !CoupledPot[r]->Saving()) CoupledPot[r]->writeData();
            S << "Coupled potential: " << (molecule != 0 ? CoupledPot[r]->getRelativePath(MolFile) : CoupledPot[r]->getFileName(true)) << '\n';
        }
    }
    Saved();
    if (writeFitData && fitData != 0)
    {
        if (Filename.isEmpty()) Filename = getFileName();
        QString FDFile = Filename.left(Filename.lastIndexOf('.')) + ".fdat", CFN = fitData->getFileName();
        QFile FDF(FDFile);
        if (!FDF.exists() || CFN != FDFile || !fitData->isSaved())
        {
            fitData->writeData(FDFile);
            if (CFN != FDFile) fitData->setFileName(CFN);
        }
        if (writePotFitTraceTab && !PotFitTraceTab.isEmpty())
        {
            QFile File(Filename.left(Filename.lastIndexOf('.')) + "_FTT.dat");
            File.open(QIODevice::WriteOnly);
            File.write((PotFitTraceTab.join("\n") + "\n").toLatin1());
            writePotFitTraceTab = false;
        }
    }
    m_saving = false;
    return true;
}

bool Potential::writePotData(QString Filename)
{
    if (!Worker->isAnaPotAv()) return false;
    double b, De, iCoeff, iOff, iExp, Ra, Ri, Rm, Tm, *LRCoeff, *PotCoeff, EI_A, EIalpha, EI_gamma;
    int nPotCoeff, nLRCoeff, *pLRCoeff, pCC = 0, n, intIExp, RIso1, RIso2;
    double Diff, D;
    if (Filename.isEmpty()) 
        if ((Filename = QFileDialog::getSaveFileName(NULL, tr("Select file name"), tr(""), tr(""))).isEmpty())
            return false;
    QFile File(Filename);
    //if (!write(&File)) return false;
    File.open(QIODevice::WriteOnly);
    QTextStream S(&File);
    QString Buffer;
    AnaPot *aPot = Worker->getAnaPotForReading();
    S << getSource() << "\n";
    aPot->getPotCoeff(b, Rm, Tm, De, nPotCoeff, PotCoeff);
    aPot->getInnerWall(Ri, iOff, iExp, iCoeff);
    aPot->getLRCoeffForReading(nLRCoeff, pLRCoeff, LRCoeff, &Ra, true);
    aPot->getExchangeCoeff(EI_A, EI_gamma, EIalpha);
    S << ("               " + QString::number(Tm, 'f', 8)).right(15);
    S << ("               " + QString::number(Rm, 'f', 8)).right(15) << "\n";
    S << ("     " + QString::number(pCC)).right(5);
    S << ("               " + QString::number(De, 'f', 8)).right(15);
    intIExp = int(iExp);
    S << ("     " + QString::number(iExp == double(intIExp) ? intIExp : -intIExp)).right(5) << "    0\n";
    S << ("     " + QString::number(nPotCoeff)).right(5);
    S << ("               " + QString::number(b, 'f', 8)).right(15);
    S << " 0.90000000D+02\n";
    for (n=0; n < nPotCoeff; n++)
    {
        S << ("     " + QString::number(n+1)).right(5);
        S << ("                         " 
                + QString::number(PotCoeff[n], 'E', 18)).right(25).replace('E', 'D') << "\n";
    }
    if (Worker->isAdCorrA())
    {
        int NAdCorr, TAdCorr, PAdCorr;
        double *adCorr, RIso1AdCorr, RIso2AdCorr, adRm, adb;
        Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, adRm, adb, 0, true);
        if (molecule != 0)
        {
            Atom *atom = molecule->getAtom1();
            for (n=1, RIso1 = atom->getnNuc(0), Diff = fabs(atom->getIsoMass(0) - RIso1AdCorr); n < atom->getnIso(); n++)
                if ((D = fabs(atom->getIsoMass(n) - RIso1AdCorr)) < Diff)
            {
                Diff = D;
                RIso1 = atom->getnNuc(n);
            }
            atom = molecule->getAtom2();
            for (n=1, RIso2 = atom->getnNuc(0), Diff = fabs(atom->getIsoMass(0) - RIso2AdCorr); n < atom->getnIso(); n++)
                if ((D = fabs(atom->getIsoMass(n) - RIso2AdCorr)) < Diff)
            {
                Diff = D;
                RIso2 = atom->getnNuc(n);
            }
        }
        else
        {
            RIso1 = lround(RIso1AdCorr);
            RIso2 = lround(RIso2AdCorr);
        }
        S << ("    " + QString::number(NAdCorr)).right(5) << ("    " + QString::number(RIso1)).right(5);
        S << ("    " + QString::number(RIso2)).right(5) << ("    " + QString::number(PAdCorr)).right(5);
        S << ("               " + QString::number(RIso1AdCorr, 'f', 9)).right(15);
        S << ("               " + QString::number(RIso2AdCorr, 'f', 9)).right(15);
        S << ("    " + QString::number(TAdCorr)).right(5) << '\n';
        for (n=0; n < NAdCorr; n++)
            S << ("    " + QString::number(n)).right(5)
              << ("                         " + QString::number(adCorr[n], 'e', 18)).right(25).replace('e', 'D') << '\n';
    }
    else S << "    0\n";
    S << "    0\n";
    S << ("     " + QString::number(nLRCoeff)).right(5) << "\n";
    for (n=0; n < nLRCoeff; n++)
    {
        S << ("     " + QString::number(pLRCoeff[n])).right(5);
        S << ("                         "
                + QString::number(LRCoeff[n], 'E', 18)).right(25).replace('E', 'D') << "\n";
    }
    if (EI_A == 0.0) EI_A = 1e-30;
    S << ("              " + QString::number(EI_A, 'e', 8)).right(15).replace('e', 'D');
    S << ("              " + QString::number(EIalpha, 'e', 8)).right(15).replace('e', 'D');
    S << ("              " + QString::number(EI_gamma, 'e', 8)).right(15).replace('e', 'D');
    S << "    0    0\n";
    S << (Buffer = ("          " + QString::number(Ri, 'f', 6)).right(10));
    S << ("          " + QString::number(Ra, 'f', 6)).right(10);
    S << "     rmin used:" << Buffer << "\n";
    S << ("                           " + QString::number(iOff, 'E', 20)).right(27).replace('E', 'D');
    S << ("                           " + 
            QString::number(iCoeff, 'E', 20)).right(27).replace('E', 'D');
    S << ("                           " + QString::number(iExp, 'E', 20)).right(27).replace('E', 'D');
    S << "\n   0.000000000000000000D+00   0.000000000000000000D+00\n";
    Worker->unlock();
    return true;
}

void Potential::writePoints(QString FileName, double RMin, double RMax, int NPoints, int NDigits, bool FirstRow, bool AdCorrE, int Iso,
                            QString MCSDir)
{
    int n, f;
    double r, h = (RMax - RMin) / double(NPoints - 1), *Dat = getPoints(RMin, RMax, NPoints);
    IsoTab *IsoT = Worker->getIsoT();
    double *adC = (AdCorrE ? getAdCorrEnergy(RMin, h, NPoints, IsoT->mIso1[Iso], IsoT->mIso2[Iso]) : 0);
    double M, E, *Min, *Max, *AMin, *AMax;
    bool MM = !MCSDir.isEmpty();
    QFile Datei(FileName);
    Datei.open(QIODevice::WriteOnly);
    QTextStream S(&Datei);
    if (MM)
    {
        Min = new double[NPoints]; Max = new double[NPoints];
        for (n=0; n < NPoints; n++) Min[n] = Max[n] = Dat[n];
        if (AdCorrE)
        {
            AMin = new double[NPoints]; AMax = new double[NPoints];
            for (n=0; n < NPoints; n++) AMin[n] = AMax[n] = adC[n];
        }
        else AMin = AMax = 0;
        QDir MDir(MCSDir);
        QFileInfoList FList = MDir.entryInfoList(QDir::Dirs | QDir::Files | QDir::NoDotAndDotDot);
        QFile MFile;
        Potential *TPot = new Potential();
        double *MDat, *MadC;
        for (f=0; f < FList.count(); f++)
        {
            if (FList[f].isDir()) 
            {
                MDir.setPath(FList[f].filePath());
                FList << MDir.entryInfoList(QDir::Dirs | QDir::Files | QDir::NoDotAndDotDot);
            }
            else if (FList[f].fileName().right(4) == ".pot")
            {
                TPot->readData(FList[f].filePath());
                MDat = TPot->getPoints(RMin, RMax, NPoints);
                for (n=0; n < NPoints; n++)
                {
                    if (MDat[n] < Min[n]) Min[n] = MDat[n];
                    else if (MDat[n] > Max[n]) Max[n] = MDat[n];
                }
                delete[] MDat;
                if (AdCorrE)
                {
                    MadC = TPot->getAdCorrEnergy(RMin, h, NPoints, IsoT->mIso1[Iso], IsoT->mIso2[Iso]);
                    for (n=0; n < NPoints; n++)
                    {
                        if (MadC[n] < AMin[n]) AMin[n] = MadC[n];
                        else if (MadC[n] > AMax[n]) AMax[n] = MadC[n];
                    }
                    delete[] MadC;
                }
            }
        }
        delete TPot;
    }
    else Min = Max = AMin = AMax = 0;
    if (FirstRow)
    {
        getMinimum(M, E);
        S << QString::number(NPoints);
        S << ("           " + QString::number(M, 'f', 6)).right(12) << "\n";
    }
    for (r = RMin, n=0; n < NPoints; n++, r+=h) 
    {
        S << ("            " + QString::number(r, 'f', 6)).right(12);
        S << ("            " + QString::number(Dat[n], 'f', NDigits)).right(NDigits + 8);
        if (MM)
        {
            S << ("            " + QString::number(Min[n], 'f', NDigits)).right(NDigits + 8);
            S << ("            " + QString::number(Max[n], 'f', NDigits)).right(NDigits + 8);
        }
        if (AdCorrE) 
        {
            S << ("            " + QString::number(adC[n], 'f', NDigits)).right(NDigits + 8);
            if (MM)
            {
                S << ("            " + QString::number(AMin[n], 'f', NDigits)).right(NDigits + 8);
                S << ("            " + QString::number(AMax[n], 'f', NDigits)).right(NDigits + 8);
            }
        }
        S << '\n';
    }
    delete[] Dat;
    if (MM)
    {
        delete[] Min;
        delete[] Max;
        if (AdCorrE)
        {
            delete[] AMin;
            delete[] AMax;
        }
    }
    if (AdCorrE) delete[] adC;
}

bool Potential::readPointPot(QString FileName)
{
    //printf("Potential::readPointPot\n");
    int n=0, N=1000, c, m=0;
    Tab->setRowCount(N);
    QFile File(FileName);
    if (!read(&File)) return false;
    QTextStream S(&File);
    QString Buffer;
    QStringList L;
    mWasMoving = false;
    Tab->blockSignals(true);
    while (!S.atEnd())
    {
        Buffer = S.readLine();
        L = Buffer.split(" ", QString::SkipEmptyParts);
        c = L.count();
        if (c!=2) 
        {
            if (n>0) if (m++ == 3) return false;
            continue;
        }
        if (L[0].toDouble() == 0.0 || (n==0 && L[1].toDouble() == 0.0)) continue;
        if (n == N) Tab->setRowCount(N+=1000);
        Tab->setItem(n, 0, new QTableWidgetItem(L[0]));
        Tab->setItem(n++, 1, new QTableWidgetItem(L[1]));
    }
    if (n < 3) return false;
    Tab->setRowCount(N = n);
    SplinePoint *points = new SplinePoint[N];
    for (n=0; n<N; n++)
    {
        points[n].x = Tab->item(n, 0)->text().toDouble();
        points[n].y = Tab->item(n, 1)->text().toDouble();
        if (n>0 ? points[n-1].x >= points[n].x : false) 
        {
            printf("read PointPot error: x[%d]=%f, x[%d]=%f\n", n-1, points[n-1].x, n, points[n].x);
            delete[] points;
            return false;
        }
    }
    if (N>2 ? points[N-1].x - points[N-2].x > 10.0 * (points[N-2].x - points[N-3].x) : false)
    {
        Worker->setUinf(points[N-1].y);
        N--;
        Tab->item(N-1, 0)->setText("U_inf");
    }
    if (!Worker->setSplinePot(N, points, 0, 0, 0, 0.0, 0.0, 0.0))
    {
        delete[] points;
        return false;
    }
    Tab->setColumnCount(3);
    Tab->setHorizontalHeaderLabels(QStringList() << "R [A]" << "Energy [cm-1]" << "d^2E/dR^2");
    Tab->blockSignals(false);
    //printf("Vor calcyss\n");
    calcyss(false);
    //printf("Ende readPointPot\n");
    return true;
}

void Potential::calcyss(const bool movingPoints)
{
    Worker->calcyss(movingPoints);
    UpdateTab();
    Changed();
    //printf("Potential::calcyss\n");
    /*int numSplinePoints, n, NLRC, *pLRC;
    SplinePoint *points;
    double *LRC, iA, iO, iExp;
    Worker->getSplinePotForWriting(numSplinePoints, points, NLRC, pLRC, LRC, iA, iO, iExp);
    if (numSplinePoints == 0) return;
    int N = numSplinePoints, st = N - 3;
    double B[N-3], F, Uinf = Worker->getUinf();
    points[0].yss = 0.0;
    if (N<3) points[N-1].yss = 0.0;
    else
    {
        if (2.0 * (points[N-2].x - points[N-3].x) < points[N-1].x - points[N-2].x 
                  && NLRC == 0 && Uinf == 0.0)
        {
            Uinf = points[N-1].x;
            N--;
            st--;
            //SLIC6 = true;
            //st = N - 4;
            //double R2 = points[N-2].x;
            //double R26 = pow(R2, 6.0), R16 = pow(points[N-1].x, 6.0);
            //points[N-1].yss = (points[N-2].y - points[N-1].y) * R26 * R16 / (R16 - R26);
            //points[N-2].yss = -42.0 * points[N-1].yss / (R26 * R2 * R2);*/
                    //0.0;
    /*    }
        points[N-1].yss = 0.0;
        if (st==0)
        { 
            points[1].yss = 3.0 / (points[2].x - points[0].x) 
                           * ((points[2].y - points[1].y) / (points[2].x - points[1].x)
                             -(points[1].y - points[0].y) / (points[1].x - points[0].x)
                             -(points[2].x - points[1].x) / 6.0 * points[2].yss);
            return;
        }
        else
        {
            F = 0.5 * (points[2].x - points[1].x) / (points[2].x - points[0].x);
            B[0] = points[3].x - points[1].x - F * 0.5 * (points[2].x - points[1].x);
            points[1].yss = 3.0 * ((points[2].y - points[1].y) / (points[2].x - points[1].x)
                                  -(points[1].y - points[0].y) / (points[1].x - points[0].x));
            points[2].yss = 3.0 * ((points[3].y - points[2].y) / (points[3].x - points[2].x)
                                  -(points[2].y - points[1].y) / (points[2].x - points[1].x))
                          - F * points[1].yss;
            for (n=1; n < st; n++)
            {
                F = 0.5 * (points[n+2].x - points[n+1].x) / B[n-1];
                B[n] = points[n+3].x - points[n+1].x - F * 0.5 * (points[n+2].x - points[n+1].x);
                points[n+2].yss = 3.0 * ((points[n+3].y - points[n+2].y) 
                                / (points[n+3].x - points[n+2].x)
                                -(points[n+2].y - points[n+1].y) / (points[n+2].x - points[n+1].x))
                                - F * points[n+1].yss;
            }
            for (n = st - 1; n >= 0; n--) 
                points[n+2].yss = 
                    (points[n+2].yss - 0.5 * (points[n+3].x - points[n+2].x) * points[n+3].yss) / B[n];
            points[1].yss = (points[1].yss - 0.5 * (points[2].x - points[1].x) * points[2].yss) 
                          / (points[2].x - points[0].x);
        }
    }
    Tab->blockSignals(true);
    for (n=0; n<N; n++) 
        Tab->setItem(n, 2, new QTableWidgetItem(QString::number(points[n].yss, 'f', 7)));
    Tab->blockSignals(false);
    Worker->recalcSCLC();
    if (N != numSplinePoints) Worker->setSplinePot(N, points, NLRC, pLRC, LRC, iA, iO, iExp);
    //calcLMatrix();
    //calcS(NumWFPoints, points[0].x, points[N-2].x);*/
}

void Potential::cdConnectLR(int p1)
{
    if (Fit->isRunning()) return;
    Worker->cdConnectLR(p1);
    UpdateTab();
    Changed();
}

void Potential::cdConnectSR()
{
    if (Fit->isRunning()) return;
    Worker->cdConnectSR(mWasMoving);
    UpdateTab();
    Changed();
}

void Potential::setPointData(int N, double *R, double *E)
{
    if (Fit->isRunning()) return;
    //printf("Potential::setPointData\n");
    //printf("R[0][1]=%f\n", R[1]);
    int n;
    Tab->setHorizontalHeaderLabels(QStringList() << "R [A]" << "E [cm^-1]" << "d^2E/dR^2");
    SplinePoint *points = new SplinePoint[N];
    Tab->setColumnCount(3);
    Tab->setRowCount(N);
    Tab->blockSignals(true);
    for (n=0; n<N; n++)
    {
        Tab->setItem(n, 0, new QTableWidgetItem(QString::number(points[n].x = R[n], 'f', 6)));
        Tab->setItem(n, 1, new QTableWidgetItem(QString::number(points[n].y = E[n], 'f', 6)));
        points[n].variable = true;
    }
    Tab->blockSignals(false);
    Worker->setSplinePot(N, points, 0, 0, 0, 0.0, 0.0, 0.0);
    calcyss(false);
    Changed();
}

void Potential::setMTTCoefficients(double B, double alpha1, int NA, double* A, double alpha, double beta, double gamma, int NLRC,
                                   int* pLRC, double* LRC, double Uinf, int NC, double *C, double Calpha, double Cbeta, 
                                   int NCLRC, int *pCLRC, double *CLRC, double CInf)
{
    if (Fit->isRunning()) return;
    MTTPot *mTTPot = new MTTPot(Fit);
    mTTPot->setCoefficients(B, alpha1, NA, A, alpha, beta, gamma, NLRC, pLRC, LRC, Uinf, NC, C, Calpha, Cbeta, NCLRC, pCLRC, 
                            CLRC, CInf);
    setWorker(mTTPot);
    UpdateTab();
    Changed();
}

bool Potential::readCoupledPotData(QString Filename)
{
    double b, De = 0.0, Ra, Ri, Rm, Tm, *LRCoeff = 0, *PotCoeff, iExp, q_e_si = 0.0, q_e_pi = 0.0, q_f_pi = 0.0, ASO_asym = 0.0, *ASO = 0, alpha = 0.0, xi_eps = 0.0, xi_Rx_s = 0.0;
    double beta = 0.0, ASO_scaling = 0.0, xi_scaling = 0.0, *xi = 0, gamma = 0.0, B_xi = 0.0, alpha_xi = 0.0, RFakt = a0_Angstrom, Rc = 0.0, Rs = 0.0, epsilon = 0.0;
    int nPotCoeff, nLRCoeff, *pLRCoeff = 0, pCC, n, m, p, N;
    QFile File(Filename);
    if (!read(&File)) return false;
    QTextStream St(&File);
    QString Buffer;
    Buffer = St.readLine();
    N = Buffer.left(5).toInt() + 2;
    if (N <= 0 || N > 5) return false;
    Potential *Pot[N];
    AnaPot *aPot;
    MTTPot *couplingFuncs;
    mWasMoving = false;
    for (p=0; p<N-2; p++)
    {
        if (p>0) 
        {
            Pot[p] = MW->CreatePotential();
            if (Pot[p] == 0) return false;
        }
        else Pot[p] = this;
        aPot = new AnaPot(Pot[p]->Fit);
        Pot[p]->setName((getFName().split('.'))[0] + QString::number(p));
        Pot[p]->setSource(Filename);
        Buffer = St.readLine().replace('D', 'e');
        Tm = Buffer.left(15).toDouble();
        Rm = Buffer.mid(15, 15).toDouble();
        Buffer = St.readLine();
        iExp = Buffer.left(5).toInt();
        De = Buffer.mid(5, 15).toDouble();
        pCC = Buffer.mid(20, 5).toInt();
        Buffer = St.readLine().replace('D', 'e');
        nPotCoeff = Buffer.left(5).toInt();
        b = Buffer.mid(5, 15).toDouble();
        PotCoeff = new double[nPotCoeff];
        for (n=0; n < nPotCoeff; n++)
        {
            Buffer = St.readLine().replace('D', 'e');
            m = Buffer.left(5).toInt();
            if (m > nPotCoeff) return false;
            PotCoeff[m-1] = Buffer.mid(5, 25).toDouble();
        }
        Buffer = St.readLine();
        nLRCoeff = Buffer.toInt();
        pLRCoeff = new int[nLRCoeff + (pCC > 0 ? 2 : 1)];
        LRCoeff = new double[nLRCoeff + (pCC > 0 ? 2 : 1)];
        for (n=0; n < nLRCoeff; n++)
        {
            Buffer = St.readLine().replace('D', 'e');
            pLRCoeff[n] = Buffer.left(5).toInt();
            LRCoeff[n] = Buffer.mid(5, 25).toDouble();
        }
        Buffer = St.readLine().replace('D', 'e');
        if (!Buffer.contains("rmin"))
        {
            int NAdCorr = Buffer.left(5).toInt();
            if (NAdCorr > 0)
            {

            }
            Buffer = St.readLine().replace('D', 'e');
        }
        Ri = Buffer.left(10).toDouble();
        Ra = Buffer.mid(10, 10).toDouble();
        aPot->setPotCoeff(b, Rm, Tm, De, nPotCoeff, PotCoeff);
        LRCoeff[n] = Buffer.mid(46, 15).toDouble();
        if (LRCoeff[n] != 0.0) nLRCoeff++;
        if (pCC < 0) pLRCoeff[n] = -pCC;
        else
        {
            pLRCoeff[n] = pCC;
            aPot->setLRCoeff(Ra, nLRCoeff, pLRCoeff, LRCoeff);
            double E1 = aPot->getPoint(Ra, 2), E2 = aPot->getPoint(Ra, 3);
            pLRCoeff[nLRCoeff] = (pCC += 2);
            LRCoeff[nLRCoeff++] = (E2 - E1) / pow(Ra, -pCC);
        }
        aPot->setLRCoeff(Ra, nLRCoeff, pLRCoeff, LRCoeff);
        aPot->setInnerWall(Ri, 0.0, iExp, 0.0);
        aPot->cdConnectSR(mWasMoving);
        Pot[p]->setWorker(aPot);
    }
    while (!St.atEnd())
    {
        Buffer = St.readLine();
        if (Buffer.indexOf("q_e sigma") >= 0) q_e_si = (Buffer.split(':'))[1].toDouble();
        else if (Buffer.indexOf("q_e pi") >= 0) q_e_pi = (Buffer.split(':'))[1].toDouble();
        else if (Buffer.indexOf("q_f pi") >= 0) q_f_pi = (Buffer.split(':'))[1].toDouble();
        else if (Buffer.indexOf("ASO_asym") >= 0) ASO_asym = (Buffer.split(':'))[1].toDouble() * hartree_cm;
        else if (Buffer.indexOf("ASO_scaling") >= 0) 
        {
            ASO_scaling = (Buffer.split(':'))[1].toDouble();
            if (ASO != 0) for (n=0; n<2; n++) ASO[n] *= ASO_scaling;
            if (LRCoeff == 0) LRCoeff[0] *= ASO_scaling;
            Pot[p] = MW->CreatePotential();
            if (Pot[p] == 0)
            {
                if (ASO != 0)
                {
                    delete[] ASO;
                    ASO = 0;
                }
                if (LRCoeff != 0)
                {
                    delete[] LRCoeff;
                    delete[] pLRCoeff;
                    LRCoeff = 0;
                    pLRCoeff = 0;
                }
                break;
            }
            couplingFuncs = new MTTPot(Pot[p]->Fit);
            couplingFuncs->setCoefficients(0.0, 0.0, 2, ASO, alpha, beta, 0.0, 1, pLRCoeff, LRCoeff, ASO_asym);
            Pot[p++]->setWorker(couplingFuncs);
        }
        else if ((m = Buffer.indexOf("ASO_")) >= 0)
        {
            if ((n = Buffer.mid(m+4, 1).toInt()) >= 0 && n<=1)
            {
                if (ASO == 0) ASO = new double[2];
                ASO[n] = (Buffer.split(':'))[1].toDouble() * hartree_cm;
                if (n==1) ASO[1] /= RFakt;
            }
        }
        else if (Buffer.indexOf("alpha_xi") >= 0) 
        {
            alpha_xi = (Buffer.split(':'))[1].toDouble() / RFakt; 
            if (xi == 0 || LRCoeff == 0) break;
            Pot[p] = MW->CreatePotential();
            if (Pot[p] == 0)
            {
                delete[] xi;
                delete[] LRCoeff;
                delete[] pLRCoeff;
                xi = 0;
                LRCoeff = 0;
                pLRCoeff = 0;
                break;
            }
            couplingFuncs = new MTTPot(Pot[p]->Fit);
            couplingFuncs->setCoefficients(B_xi, alpha_xi, 3, xi, 0.0, beta, gamma, 1, pLRCoeff, LRCoeff, 0.0);
            Pot[p++]->setWorker(couplingFuncs);
        }
        else if (Buffer.indexOf("alpha") >= 0) alpha = (Buffer.split(':'))[1].toDouble() / RFakt;
        else if (Buffer.indexOf("beta") >= 0) beta = (Buffer.split(':'))[1].toDouble() / RFakt;
        else if (Buffer.indexOf("SO_C6") >= 0)
        {
            pLRCoeff = new int[1];
            pLRCoeff[0] = 6;
            LRCoeff = new double[1];
            LRCoeff[0] = -(Buffer.split(':'))[1].toDouble() * hartree_cm * pow(RFakt, 6);
        }
        else if (Buffer.indexOf("xi_scaling") >= 0) xi_scaling = (Buffer.split(':'))[1].toDouble();
        else if (Buffer.indexOf("gamma") >= 0) gamma = (Buffer.split(':'))[1].toDouble() / (RFakt * RFakt);
        else if (Buffer.indexOf("xi_C4") >= 0)
        {
            LRCoeff = new double[1];
            LRCoeff[0] = -(Buffer.split(':'))[1].toDouble() * hartree_cm * pow(RFakt, 4) * xi_scaling;
            pLRCoeff = new int[1];
            pLRCoeff[0] = 4;
        }
        else if (Buffer.indexOf("B_xi") >= 0) B_xi = (Buffer.split(':'))[1].toDouble() * hartree_cm * xi_scaling;
        else if (Buffer.indexOf("xi_eps") >= 0) xi_eps = (Buffer.split(':'))[1].toDouble();
        else if (Buffer.indexOf("xi_Rx_s") >= 0) xi_Rx_s = (Buffer.split(':'))[1].toDouble();
        else if ((m = Buffer.indexOf("xi_")) >= 0)
        {
            if ((n = Buffer.mid(m+3, 1).toInt()) >= 0 && n<=1)
            {
                if (xi == 0) xi = new double[2];
                xi[n] = (Buffer.split(':'))[1].toDouble() * (xi_scaling != 0.0 ? xi_scaling * (n==0 ? hartree_cm : hartree_cm / (n==1 ? RFakt : RFakt * RFakt)) : 1.0);
            }
        }
        else if (Buffer.left(15).toDouble() != 0.0 && Buffer.mid(15, 15).toDouble() != 0.0)
        {
            Rc = Buffer.left(15).toDouble();
            Rs = Buffer.mid(15, 15).toDouble();
            epsilon = Buffer.mid(30, 15).toDouble();
        }
    }
    if (xi_scaling == 0 && Rc != 0 && Rs != 0.0 && epsilon != 0.0)
    {
        m_couplingFuncs = new CouplingFuncs(MW);
        m_couplingFuncs->setData(2, xi, 0.0, Rs, Rc, epsilon, xi_Rx_s, xi_eps);
        MW->showMDIChild(m_couplingFuncs);
    }
    else m_couplingFuncs = 0;
    if ((N = p) < 2) return false;
    if (q_e_si != 0.0) Worker->setLambdaDoubling(q_e_si, 0.0);
    if (q_e_pi != 0.0 || q_f_pi != 0.0) Pot[1]->Worker->setLambdaDoubling(q_e_pi, q_f_pi);
    if (CoupledPot != 0)
    {
        delete[] CoupledPot;
        delete[] CoupledComp;
    }
    CoupledPot = new Potential*[NCoupled = N - 1];
    CoupledComp = new int[NCoupled];
    for (p=0; p<N; p++)
    {
        Pot[p]->UpdateTab();
        Pot[p]->setImported();
        Pot[p]->Tab->setHorizontalHeaderLabels(QStringList() << "coefficient" << "value" << "error" 
                    << "sig. digits");
        Pot[p]->Changed();
        if (p>0) CoupledPot[p-1] = Pot[p];
        Pot[p]->m_couplingFuncs = m_couplingFuncs;
        Pot[p]->show();
    }
    for (p=1; p < N; ++p)
    {
        Pot[p]->NCoupled = NCoupled;
        Pot[p]->CoupledPot = new Potential*[NCoupled];
        Pot[p]->CoupledComp = new int[NCoupled];
        for (n=m=0; m < N; ++m) if (p!=m) Pot[p]->CoupledPot[n++] = Pot[m];
    }
    setImported();
    Changed();
    return true;
}

bool Potential::readPotData(QString FileName)
{
    if (Fit->isRunning()) return false;
    double De = 0.0, iCoeff, iOff, Ra, Ri, Tm, *LRCoeff, *PotCoeff, iExp, ExcInt_A, ExcInt_Alpha, ExcInt_gamma, AdCorr_Rm, AdCorr_b;
    double RIso1AdCorr, RIso2AdCorr, *adCorr;
    bool *adCorrFree;
    int nPotCoeff, nLRCoeff, *pLRCoeff, pCC, n, m, NAdCorr, PAdCorr, TAdCorr;
    QFile File(FileName);
    if (!read(&File)) return false;
    QTextStream St(&File);
    QString Buffer;
    AnaPot *aPot = new AnaPot(Fit);
    mWasMoving = false;
    Buffer = St.readLine();
    if (Buffer.left(5).toInt() > 0) return false;
    setName((getFName().split('.'))[0]);
    setSource(Buffer);
    Buffer = St.readLine().replace('D', 'e');
    Tm = Buffer.left(17).toDouble();
    AdCorr_Rm = Buffer.right(Buffer.length() - 17).toDouble();
    Buffer = St.readLine();
    pCC = Buffer.left(5).toInt();
    De = Buffer.mid(5, 15).toDouble();
    if ((iExp = Buffer.mid(20, 5).toInt()) > 50.0) iExp -= 50.0;
    Buffer = St.readLine().replace('D', 'e');
    nPotCoeff = Buffer.left(5).toInt();
    AdCorr_b = Buffer.mid(5, 15).toDouble();
    if (nPotCoeff == 0) return false;
    PotCoeff = new double[nPotCoeff];
    for (n=0; n < nPotCoeff; n++)
    {
        Buffer = St.readLine().replace('D', 'e');
        m = Buffer.left(5).toInt();
        PotCoeff[m-1] = Buffer.mid(5, 25).toDouble();
    }
    Buffer = St.readLine();
    NAdCorr = Buffer.left(5).toInt();
    PAdCorr = Buffer.mid(15, 5).toInt();
    RIso1AdCorr = Buffer.mid(20, 15).toDouble();
    RIso2AdCorr = Buffer.mid(35, 15).toDouble();
    TAdCorr = Buffer.mid(50, 5).toInt();
    if (NAdCorr > 0)
    {
        adCorr = new double[NAdCorr];
        adCorrFree = new bool[NAdCorr];
        for (n=0; n < NAdCorr; n++) 
        {
            Buffer = St.readLine().replace('D', 'e');
            adCorr[n] = Buffer.right(Buffer.length() - 5).toDouble();
            adCorrFree[n] = true;
        }
    }
    else
    {
        adCorr = 0;
        adCorrFree = 0;
    }
    St.readLine();
    Buffer = St.readLine();
    nLRCoeff = Buffer.toInt();
    if (nLRCoeff == 0) return false;
    pLRCoeff = new int[nLRCoeff + 2];
    LRCoeff = new double[nLRCoeff + 2];
    for (n=0; n < nLRCoeff; n++)
    {
        Buffer = St.readLine().replace('D', 'e');
        pLRCoeff[n] = Buffer.left(5).toInt();
        LRCoeff[n] = Buffer.mid(5, 25).toDouble();
    }
    Buffer = St.readLine().replace('D', 'e');
    ExcInt_A = Buffer.left(15).toDouble();
    ExcInt_Alpha = Buffer.mid(30, 15).toDouble();
    ExcInt_gamma = Buffer.mid(15, 15).toDouble();
    Buffer = St.readLine().replace('D', 'e');
    Ri = Buffer.left(10).toDouble();
    Ra = Buffer.mid(10, 10).toDouble();
    if (Ri < 0.0 || Ri >= Ra) return false;
    Buffer = St.readLine().replace('D', 'e');
    iOff = Buffer.left(27).toDouble();
    iCoeff = Buffer.mid(27, 27).toDouble();
    if (iExp < 0.0) iExp = Buffer.mid(54, 27).toDouble();
    Buffer = St.readLine().replace('D', 'e');
    pLRCoeff[n] = pCC;
    pLRCoeff[n+1] = pCC + 2;
    LRCoeff[n] = Buffer.left(27).toDouble();
    LRCoeff[n+1] = Buffer.mid(27, 27).toDouble();
    //printf("readPotData: iExp = %f\n", iExp);
    if (LRCoeff[n] != 0.0 && LRCoeff[n+1] != 0.0) nLRCoeff += 2;
    //printf("Buffer=%s, iOff=%f, iCoeff=%f\n", Buffer.ascii(), iOff, iCoeff);
    //printf("Left=%s, Mid=%s\n", Buffer.left(27).ascii(), Buffer.mid(27,27).ascii());
    aPot->setInnerWall(Ri, iOff, iExp, iCoeff);
    aPot->setPotCoeff(AdCorr_b, AdCorr_Rm, Tm, De, nPotCoeff, PotCoeff);
    aPot->setLRCoeff(Ra, nLRCoeff, pLRCoeff, LRCoeff);
    aPot->setExchangeCoeff(ExcInt_A, ExcInt_Alpha, ExcInt_gamma);
    //aPot->sdConnectLR(pCC);
    UpdateTab();
    setWorker(aPot);
    Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
    Worker->setLambdaDoubling(0.0, 0.0);
    
    /*double Rmin = 3.3, Rmax = 100.0;
    int v, Np = 100000;
    double *R = aPot->getPointRep(Rmin, Rmax, Np), U[Np];
    if (R==NULL) 
    {
        printf("Fehler: U=NULL!\n");
        return false;
    }
    double h = (Rmax - Rmin) / (Np - 1);//, R;//[Np];
    CRotPot(19, 0, R, U, Rmin, h, Np);
    
    QFile DDatei("DebugPot.dat");
    DDatei.open(QIODevice::WriteOnly);
    QTextStream S1(&DDatei);
    for (n=1, R[0]=Rmin; n<Np; n++) //R[n] = R[n-1] + h;
        S1 << QString::number(R+=h, 'g', 10) << "    " << QString::number(U[n], 'g', 10) << "\n";
    for (n=0, R=Rmin; R <= 10.0; n++) 
    {
        if (R > 9.5) printf("R=%f, V=%f\n", R += h, U[n]);
        R+=h;
    }*/
    //double DE;//, *WF = new double[Np];
    //Numerov(9, h, -500.0, DE, v, Np, U, WF);
    //printf("v=%d\n", v);
    /*double **WF = Create(Maxv, Np), *E = new double[Maxv];
    printf("Vor NumerovCooley\n");
    NumerovCooley(0, h, n, E, WF, Np, U);
    for (m=0; m<n; m++) printf("v=%d, E=%f\n", m, E[m]);
    for (n=0; R[n] < Ri; n++);
    for (m=n; R[m] < Ra; m++);
    m-=n;
    aPot->calcPotCoeff(m, R+n, U+n, nPotCoeff, Rm, b);*/
    /*delete[] E;
    //delete[] U;
    //delete[] WF;
    delete[] R;
    Destroy(WF, Maxv);*/
    setImported();
    Tab->setHorizontalHeaderLabels(QStringList() << "coefficient" << "value" << "error" 
            << "sig. digits");
    Changed();
    return true;
}

void Potential::UpdateTab(PotentialData *D)
{
    //printf("Potential::UpdateTab\n");
    LRCTabOffs = adCorrTabOffs = ExcIntATabPos = -1;
    if (D == 0) D = Worker->getPotentialData();
    Tab->blockSignals(true);
    if (Worker->isAnaPotAv()) 
    {
        int n, r;
        QString Buffer;
        Tab->setRowCount((D->q_f != 0.0 ? 12 : 11) + (D->q_e != 0.0 ? D->nPotCoeff + 1 : D->nPotCoeff) + (D->NAdCorr > 0 ? D->NAdCorr + 4 : 0)
                            + (D->A_ex == 0.0 ? D->nLRCoeff : D->nLRCoeff + 3) + (D->NSpinRGamma > 0 ? D->NSpinRGamma + 2 : 0));
        for (n=1, r=0; n <= D->nPotCoeff; n++, r++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem((D->CFree[r] ? QIcon() : *FixPix), 
                                                     "a" + QString::number(n) + " [cm^-1]"));
            Tab->setItem(r, 1, new QTableWidgetItem((D->CFree[r] ? QIcon() : *FixPix), 
                                                    QString::number(D->PotCoeff[r], 'E', 16)));
        }
        Tab->setItem(r, 0, new QTableWidgetItem("b"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->b, 'f', 4)));
        Tab->setItem(r, 0, new QTableWidgetItem("Rm [A]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Rm, 'f', 10)));
        Tab->setItem(r, 0, new QTableWidgetItem("Tm [cm^-1]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Tm, 'f', 10)));
        Tab->setItem(r, 0, new QTableWidgetItem("Ri [A]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Ri, 'f', 10)));
        Tab->setItem(r, 0, new QTableWidgetItem("n"));
        Tab->setItem(r++, 1, new QTableWidgetItem(Buffer = QString::number(D->iExp, 'g', 10)));
        //printf("UpdateTab: iExp=%f\n", iExp);
        Tab->setItem(r, 0, new QTableWidgetItem("A [cm^-1]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->iOff, 'f', 10)));
        Tab->setItem(r, 0, new QTableWidgetItem("B [cm^-1 A^" + Buffer + "]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->iCoeff, 'f', 16)));
        Tab->setItem(r, 0, new QTableWidgetItem("Ra [A]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Ra, 'f', 10)));
        LRCTabOffs = r;
        for (n=0; n < D->nLRCoeff; n++)
        {
            Buffer = QString::number(D->pLRCoeff[n]);
            if (n < D->nLRCoeff - 1 ? D->LRCFree[n] : true)
            {
                Tab->setItem(r, 0, new QTableWidgetItem("C" + Buffer + " [cm^-1 A^" + Buffer + "]"));
                Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->LRCoeff[n], 'E', 16)));
            }
            else
            {
                Tab->setItem(r, 0, new QTableWidgetItem(*FixPix, "C" + Buffer + " [cm^-1 A^" + Buffer + "]"));
                Tab->setItem(r++, 1, new QTableWidgetItem(*FixPix, QString::number(D->LRCoeff[n], 'E', 16)));
            }
        }
        Tab->setItem(r, 0, new QTableWidgetItem("U infinity [cm^-1]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->De, 'f', 10)));
        if (D->A_ex != 0.0)
        {
            Tab->setItem(ExcIntATabPos = r, 0, new QTableWidgetItem((D->EIFree ? QIcon() : *FixPix),
                         "A_ex [cm^-1 A^" + QString::number(D->gamma, 'f', 5) + "]"));
            Tab->setItem(r++, 1, new QTableWidgetItem((D->EIFree ? QIcon() : *FixPix), QString::number(D->A_ex, 'E', 8)));
            Tab->setItem(r, 0, new QTableWidgetItem("gamma"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->gamma, 'f', 5)));
            Tab->setItem(r, 0, new QTableWidgetItem("beta [A^-1]"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->beta, 'f', 5)));
        }
        if (D->q_e != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("q_e"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->q_e, 'E', 16)));
        }
        if (D->q_f != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("q_f"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->q_f, 'E', 16)));
        }
        if (D->NAdCorr > 0) 
        {
            adCorrTabOffs = r;
            for (n=0; n < D->NAdCorr; n++)
            {
                Tab->setItem(r, 0, new QTableWidgetItem((D->adCorrFree[n] ? QIcon() : *FixPix), "AdCorr" + QString::number(n)));
                Tab->setItem(r++, 1, new QTableWidgetItem((D->adCorrFree[n] ? QIcon() : *FixPix), QString::number(D->adCorr[n], 'E', 16)));
            }
            Tab->setItem(r, 0, new QTableWidgetItem("Type AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->TAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("Pow AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->PAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso1AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso1AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso2AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso2AdCorr, 'f', 9)));
        }
        if (D->NSpinRGamma > 0)
        {
            for (n=0; n < D->NSpinRGamma; n++)
            {
                Tab->setItem(r, 0, new QTableWidgetItem("SpinRGamma" + QString::number(n)));
                Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->SpinRGamma[n], 'E', 16)));
            }
            Tab->setItem(r, 0, new QTableWidgetItem("SpinRR1"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->SpinRR1, 'E', 10)));
            Tab->setItem(r, 0, new QTableWidgetItem("SpinRR2"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->SpinRR2, 'E', 10)));
        }
        Tab->setItem(r, 0, new QTableWidgetItem("R min [A]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(rmin, 'f', 10)));
        Tab->setItem(r, 0, new QTableWidgetItem("R max [A]"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(rmax, 'f', 10)));
    }
    else if (Worker->isTTPot())
    {
        int r, n;
        Tab->setRowCount(D->N + D->NA + (D->NAdCorr > 0 ? D->NAdCorr + 8 : 2));
        for (n=0; n < D->NA; n++)
        {
            Tab->setItem(n, 0, new QTableWidgetItem("A_" + QString::number(n)));
            Tab->setItem(n, 1, new QTableWidgetItem(QString::number(D->A[n], 'f', 10)));
        }
        Tab->setItem(n, 0, new QTableWidgetItem("b"));
        Tab->setItem(n, 1, new QTableWidgetItem(QString::number(D->b, 'f', 10)));
        LRCTabOffs = r = D->NA + 1;
        for (n=0; n < D->N; n++, r++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("C_" + QString::number(2*(n+3))));
            Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->C[n], 'E', 16)));
        }
        Tab->setItem(r, 0, new QTableWidgetItem("U_inf"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Uinf, 'f', 10)));
        if (D->NAdCorr > 0)
        {
            adCorrTabOffs = r;
            for (n=0; n < D->NAdCorr; n++, r++)
            {
                Tab->setItem(r, 0, new QTableWidgetItem("AdCorr" + QString::number(n)));
                Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->adCorr[n], 'E', 16)));
            }
            Tab->setItem(r, 0, new QTableWidgetItem("Type AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->TAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("Pow AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->PAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso1AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso1AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso2AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso2AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr b"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->AdCorr_b, 'f', 4)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr Rm"));
            Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->AdCorr_Rm, 'f', 10)));
        }
    }
    else if (Worker->isLPot())
    {
        int n, r;
        Tab->setRowCount(D->NCi + D->Nb + (D->NAdCorr > 0 ? D->NAdCorr + 12 : 6));
        if (Tab->item(0, 0) != 0) Tab->item(0, 0)->setText("p");
        else Tab->setItem(0, 0, new QTableWidgetItem("p"));
        if (Tab->item(0, 1) != 0) Tab->item(0, 1)->setText(QString::number(D->p));
        else Tab->setItem(0, 1, new QTableWidgetItem(QString::number(D->p)));
        if (Tab->item(1, 0) != 0) Tab->item(1, 0)->setText("q");
        else Tab->setItem(1, 0, new QTableWidgetItem("q"));
        if (Tab->item(1, 1) != 0) Tab->item(1, 1)->setText(QString::number(D->q));
        else Tab->setItem(1, 1, new QTableWidgetItem(QString::number(D->q)));
        if (Tab->item(2, 0) != 0) Tab->item(2, 0)->setText("R_ref [A]");
        else Tab->setItem(2, 0, new QTableWidgetItem("R_ref [A]"));
        if (Tab->item(2, 1) != 0) Tab->item(2, 1)->setText(QString::number(D->Rref, 'f', 2));
        else Tab->setItem(2, 1, new QTableWidgetItem(QString::number(D->Rref, 'f', 2)));
        if (Tab->item(3, 0) != 0) Tab->item(3, 0)->setText("R_e [A]");
        else Tab->setItem(3, 0, new QTableWidgetItem("R_e [A]"));
        if (Tab->item(3, 1) != 0) Tab->item(3, 1)->setText(QString::number(D->Re, 'f', 8));
        else Tab->setItem(3, 1, new QTableWidgetItem(QString::number(D->Re, 'f', 8)));
        if (Tab->item(4, 0) != 0) Tab->item(4, 0)->setText("D_e [cm^-1]");
        else Tab->setItem(4, 0, new QTableWidgetItem("D_e [cm^-1]"));
        if (Tab->item(4, 1) != 0) Tab->item(4, 1)->setText(QString::number(D->De, 'f', 4));
        else Tab->setItem(4, 1, new QTableWidgetItem(QString::number(D->De, 'f', 4)));
         if (Tab->item(5, 0) != 0) Tab->item(5, 0)->setText("U_infinity [cm^-1]");
        else Tab->setItem(5, 0, new QTableWidgetItem("U_infinity [cm^-1]"));
        if (Tab->item(5, 1) != 0) Tab->item(5, 1)->setText(QString::number(D->T, 'f', 4));
        else Tab->setItem(5, 1, new QTableWidgetItem(QString::number(D->T, 'f', 4)));
        for (n=0; n < D->Nb; n++)
        {
            if (Tab->item(n+6, 0) != 0) 
                Tab->item(n+6, 0)->setText("beta_" + QString::number(n));
            else Tab->setItem(n+6, 0, 
                              new QTableWidgetItem("beta_" + QString::number(n)));
            if (Tab->item(n+6, 1) != 0) Tab->item(n+6, 1)->setText(QString::number(D->bA[n], 'f', 8));
            else Tab->setItem(n+6, 1, new QTableWidgetItem(QString::number(D->bA[n], 'f', 8)));
        }
        LRCTabOffs = n + D->Nb + 6;
        for (n=0; n < D->NCi; n++)
        {
            if (Tab->item(n + D->Nb + 6, 0) != 0) 
                Tab->item(n + D->Nb + 6, 0)->setText("C_" + QString::number(D->pCi[n]) + " [cm^-1 A^"
                                                  + QString::number(D->pCi[n]) + "]");
            else Tab->setItem(n + D->Nb + 6, 0, new QTableWidgetItem("C_" + QString::number(D->pCi[n]) 
                                                    + " [cm^-1 A^" + QString::number(D->pCi[n]) + "]"));
            if (Tab->item(n + D->Nb + 6, 1) != 0) 
                Tab->item(n + D->Nb + 6, 1)->setText(QString::number(D->Ci[n], 'e', 10));
            else Tab->setItem(n + D->Nb + 6, 1, new QTableWidgetItem(QString::number(D->Ci[n], 'e', 10)));
        }
        if (D->NAdCorr > 0) 
        {
            adCorrTabOffs = r = n + D->Nb + 7;
            for (n=0; n < D->NAdCorr; n++, r++)
            {
                Tab->setItem(r, 0, new QTableWidgetItem("AdCorr" + QString::number(n)));
                Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->adCorr[n], 'E', 16)));
            }
            Tab->setItem(r, 0, new QTableWidgetItem("Type AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->TAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("Pow AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->PAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso1AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso1AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso2AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso2AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr b"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->AdCorr_b, 'f', 4)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr Rm"));
            Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->AdCorr_Rm, 'f', 10)));
        }
    }
    else if (Worker->isMTTPot())
    {
        int n, r=0;
        Tab->setRowCount(D->NA + 7 + D->NLRC + D->NSA + D->NSLRC);
        if (D->B != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("B"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->B, 'g', 15)));
        }
        if (D->alpha1 != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("alpha1"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->alpha1, 'f', 14)));
        }
        for (n=0; n < D->NA; n++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("A" + QString::number(n)));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->A[n], 'g', 15)));
        }
        Tab->setItem(r, 0, new QTableWidgetItem("alpha"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->alpha, 'f', 14)));
        Tab->setItem(r, 0, new QTableWidgetItem("beta"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->beta, 'f', 14)));
        if (D->gamma != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("gamma"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->gamma, 'f', 14)));
        }
        LRCTabOffs = r;
        for (n=0; n < D->NLRC; n++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("C" + QString::number(D->pLRC[n])));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->LRC[n], 'g', 15)));
        }
        Tab->setItem(r, 0, new QTableWidgetItem("U_inf"));
        Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Uinf, 'f', 4)));
        for (n=0; n < D->NSA; n++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("SA" + QString::number(n)));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Salpha, 'f', 14)));
        }
        if (D->Salpha != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("Salpha"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Salpha, 'f', 14)));
        }
        if (D->Sbeta != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("Sbeta"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Sbeta, 'f', 14)));
        }
        for (n=0; n < D->NSLRC; n++)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("SC" + QString::number(D->pSLRC[n])));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->SLRC[n], 'g', 15)));
        }
        if (D->Sinf != 0.0)
        {
            Tab->setItem(r, 0, new QTableWidgetItem("S_inf"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->Sinf, 'f', 4)));
        }
        if (D->NAdCorr > 0) 
        {
            adCorrTabOffs = r;
            for (n=0; n < D->NAdCorr; n++, r++)
            {
                Tab->setItem(r, 0, new QTableWidgetItem("AdCorr" + QString::number(n)));
                Tab->setItem(r, 1, new QTableWidgetItem(QString::number(D->adCorr[n], 'E', 16)));
            }
            Tab->setItem(r, 0, new QTableWidgetItem("Type AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->TAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("Pow AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->PAdCorr)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso1AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso1AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("RIso2AdCorr"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->RIso2AdCorr, 'f', 9)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr b"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->AdCorr_b, 'f', 4)));
            Tab->setItem(r, 0, new QTableWidgetItem("AdCorr Rm"));
            Tab->setItem(r++, 1, new QTableWidgetItem(QString::number(D->AdCorr_Rm, 'f', 10)));
        }
        Tab->setRowCount(r);
    }
    else if (Worker->isSplinePot())
    {
        int c, n;
        Tab->setRowCount(D->N + D->NLRC + (D->NAdCorr > 0 ? D->NAdCorr + 6 : 0) 
                          + (D->Ro == 0 ? (D->iExp != 6.0 ? 4 : 3) : 2));
        for (n=0; n < D->N; n++)
        {
            if (D->points[n].variable)
            {
                Tab->setItem(n, 0, new QTableWidgetItem(QString::number(D->points[n].x, 'f', 6)));
                Tab->setItem(n, 1, new QTableWidgetItem(QString::number(D->points[n].y, 'f', 6)));
                Tab->setItem(n, 2, new QTableWidgetItem(QString::number(D->points[n].yss, 'f', 6)));
            }
            else
            {
                Tab->setItem(n, 0, new QTableWidgetItem(*FixPix, QString::number(D->points[n].x, 'f', 6)));
                Tab->setItem(n, 1, new QTableWidgetItem(*FixPix, QString::number(D->points[n].y, 'f', 6)));
                Tab->setItem(n, 2, new QTableWidgetItem(*FixPix, QString::number(D->points[n].yss, 'f', 6)));
            }
        }
        if (D->Ro == 0.0)
        {
            if (D->shortRVariable)
            {
                Tab->setItem(n, 0, new QTableWidgetItem("iA"));
                Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->iA, 'e', 12)));
                Tab->setItem(n, 0, new QTableWidgetItem("iO"));
                Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->iO, 'f', 4)));
                if (D->iExp != 6.0)
                {
                    Tab->setItem(n, 0, new QTableWidgetItem("iExp"));
                    Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->iExp, 'f', 10)));
                }
            }
            else
            {
                Tab->setItem(n, 0, new QTableWidgetItem(*FixPix, "iA"));
                Tab->setItem(n++, 1, new QTableWidgetItem(*FixPix, QString::number(D->iA, 'e', 12)));
                Tab->setItem(n, 0, new QTableWidgetItem(*FixPix, "iO"));
                Tab->setItem(n++, 1, new QTableWidgetItem(*FixPix, QString::number(D->iO, 'f', 4)));
                if (D->iExp != 6.0)
                {
                    Tab->setItem(n, 0, new QTableWidgetItem(*FixPix, "iExp"));
                    Tab->setItem(n++, 1, new QTableWidgetItem(*FixPix, QString::number(D->iExp, 'f', 10)));
                }
            }
        }
        else
        {
            Tab->setItem(n, 0, new QTableWidgetItem("R_o"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->Ro, 'f', 4)));
        }
        Tab->setItem(n, 0, new QTableWidgetItem("U_inf"));
        Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->Uinf, 'f', 4)));
        LRCTabOffs = n;
        for (c=0; c < D->NLRC; c++)
        {
            Tab->setItem(n, 0, 
                         new QTableWidgetItem((D->LRCFree[c] ? QIcon() : *FixPix), 
                                              "C_" + QString::number(D->PLRC[c])));
            Tab->setItem(n++, 1, 
                         new QTableWidgetItem((D->LRCFree[c] ? QIcon() : *FixPix),
                                              QString::number(D->LRC[c], 'e', 12)));
        }
        if (D->NAdCorr > 0) 
        {
            adCorrTabOffs = n;
            for (c=0; c < D->NAdCorr; n++, c++)
            {
                Tab->setItem(n, 0, 
                             new QTableWidgetItem((D->adCorrFree[c] ? QIcon() : *FixPix),
                                                  "AdCorr" + QString::number(c)));
                Tab->setItem(n, 1, 
                             new QTableWidgetItem((D->adCorrFree[c] ? QIcon() : *FixPix),
                                                  QString::number(D->adCorr[c], 'E', 16)));
            }
            Tab->setItem(n, 0, new QTableWidgetItem("Type AdCorr"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->TAdCorr)));
            Tab->setItem(n, 0, new QTableWidgetItem("Pow AdCorr"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->PAdCorr)));
            Tab->setItem(n, 0, new QTableWidgetItem("RIso1AdCorr"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->RIso1AdCorr, 'f', 9)));
            Tab->setItem(n, 0, new QTableWidgetItem("RIso2AdCorr"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->RIso2AdCorr, 'f', 9)));
            Tab->setItem(n, 0, new QTableWidgetItem("AdCorr b"));
            Tab->setItem(n++, 1, new QTableWidgetItem(QString::number(D->AdCorr_b, 'f', 4)));
            Tab->setItem(n, 0, new QTableWidgetItem("AdCorr Rm"));
            Tab->setItem(n, 1, new QTableWidgetItem(QString::number(D->AdCorr_Rm, 'f', 10)));
        }
    }
    //printf("Ende UpdateTab\n");
    Tab->blockSignals(false);
    delete D;
}

void Potential::VaryCoefficients(bool V, int NRows, int *Rows)
{
    if (!(Worker->isAnaPotAv() || Worker->isSplinePot())) 
    {
        printf("Potential::VaryCoefficients is not implemented for the case numSplinePoints == 0 && aPot == 0 !\n");
        return;
    }
    if (Fit->isRunning()) return;
    int i, r, N;
    if (NRows == 0)
    {
        QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
        for (i=N=0; i < SR.count(); i++) for (r=SR[i].topRow(); r <= SR[i].bottomRow(); r++) N++;
        Rows = new int[N];
        for (i=N=0; i < SR.count(); i++) for (r=SR[i].topRow(); r <= SR[i].bottomRow(); r++) Rows[N++] = r;
    }
    else N = NRows;
    Tab->blockSignals(true);
    if (Worker->isSplinePot())
    {
        int numSplinePoints, NLRC, *PLRC, NAdCorr, TAdCorr, PAdCorr;
        SplinePoint *points;
        double *LRC, iA, iO, iExp, *adCorr, RIso1, RIso2, Rm, b;
        bool shortRVariable = Worker->getShortRVariable(), *LRCfree = Worker->getLRCFree(), *adCorrFree;
        Worker->getSplinePotForWriting(numSplinePoints, points, NLRC, PLRC, LRC, iA, iO, iExp);
        Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b, &adCorrFree, true);
        int NSR = (iExp == 6.0 ? 2 : 3), c;
        //bool GS = (State != 0 ? State->getStateNum() == 0 : false);
        for (r=0; r<N; r++)
        {
            if (Rows[r] < numSplinePoints) 
            {
                if ((points[Rows[r]].variable = V)) 
                {
                    Tab->item(Rows[r], 0)->setIcon(QIcon());
                    Tab->item(Rows[r], 1)->setIcon(QIcon());
                }
                else
                {
                    Tab->item(Rows[r], 0)->setIcon(*FixPix);
                    Tab->item(Rows[r], 1)->setIcon(*FixPix);
                }
            }
            else if (Rows[r] < numSplinePoints + NSR)
            {
                if ((shortRVariable = V))
                {
                    Tab->item(numSplinePoints, 0)->setIcon(QIcon());
                    Tab->item(numSplinePoints, 1)->setIcon(QIcon());
                    Tab->item(numSplinePoints + 1, 0)->setIcon(QIcon());
                    Tab->item(numSplinePoints + 1, 1)->setIcon(QIcon());
                }
                else
                {
                    Tab->item(numSplinePoints, 0)->setIcon(*FixPix);
                    Tab->item(numSplinePoints, 1)->setIcon(*FixPix);
                    Tab->item(numSplinePoints + 1, 0)->setIcon(*FixPix);
                    Tab->item(numSplinePoints + 1, 1)->setIcon(*FixPix);
                }
                //r = numSplinePoints + NSR;
            }
            else if (Rows[r] < numSplinePoints + NSR + NLRC + 1)
            {
                //if (GS)
                //{
                    if ((c = Rows[r] - numSplinePoints - NSR - 1) >= 0)
                    {
                        if ((LRCfree[c] = V))
                        {
                            Tab->item(Rows[r], 0)->setIcon(QIcon());
                            Tab->item(Rows[r], 1)->setIcon(QIcon());
                        }
                        else
                        {
                            Tab->item(Rows[r], 0)->setIcon(*FixPix);
                            Tab->item(Rows[r], 1)->setIcon(*FixPix);
                        }
                    }
                /*}
                else
                {
                    if (V) for (r = numSplinePoints + NSR; r < numSplinePoints + NSR + NLRC + 1; r++)
                    {
                        Tab->item(r, 0)->setIcon(QIcon());
                        Tab->item(r, 1)->setIcon(QIcon());
                    }
                    else for (r = numSplinePoints + 2; r < Tab->rowCount(); r++)
                    {
                        Tab->item(r, 0)->setIcon(*FixPix);
                        Tab->item(r, 1)->setIcon(*FixPix);
                    }
                    for (c=0; c < NLRC; c++) LRCfree[c] = V;
                }*/
            }
            else if ((c = Rows[r] - numSplinePoints - NSR - NLRC - 1) < NAdCorr)
            {
                if ((adCorrFree[c] = V))
                {
                    Tab->item(Rows[r], 0)->setIcon(QIcon());
                    Tab->item(Rows[r], 1)->setIcon(QIcon());
                }
                else
                {
                    Tab->item(Rows[r], 0)->setIcon(*FixPix);
                    Tab->item(Rows[r], 1)->setIcon(*FixPix);
                }
            }
        }
    }
    else
    {
        for (r=0; r<N; r++)
        {
            Tab->item(Rows[r], 0)->setIcon(V ? QIcon() : *FixPix);
            Tab->item(Rows[r], 1)->setIcon(V ? QIcon() : *FixPix);
        }
        updatePot();
    }
    if (N>0) 
    {
        if (NRows == 0) delete[] Rows;
        Worker->DestroySMap();
    }
    Tab->blockSignals(false);
}

void Potential::getFixedCoefficients(int& N, int*& Rows)
{
    int n, NR = Tab->rowCount();
    for (n=N=0; n < NR; n++) if (!Tab->item(n, 0)->icon().isNull()) N++;
    Rows = new int[N];
    for (n=N=0; n < NR; n++) if (!Tab->item(n, 0)->icon().isNull()) Rows[N++] = n;
}

void Potential::updatePot()
{
    UpdatePot(Worker->getType() != NoPotential ? unknown : likelyASplinePotential);
    Changed();
}

void Potential::UpdatePot(PreliminaryPotentialType type)
{
    if (Fit->isRunning()) return;
    mWasMoving = false;
    if (type == likelyASplinePotential)
    {
        if (testIfSplinePot())
        {
            type = splinePot;
            /*if (2.0 * (points[N-2].x - points[N-3].x) < points[N-1].x - points[N-2].x 
                         && points[N-1].yss != 0.0) 
            {
                printf("SLIC6!\n");
                SLIC6 = true;
            }*/
        }
        else
        {
            Tab->setHorizontalHeaderLabels(QStringList() 
                    << "coefficient" << "value" << "error" << "sig. digits");
            if (Tab->item(2, 0)->text().indexOf("R_ref") >= 0) type = morseLongRange;
            else if (Tab->item(0, 0)->text().indexOf("A0") >= 0) type = modTangToen;
            else type = anaPot;
        }
    }
    LRCTabOffs = adCorrTabOffs = ExcIntATabPos = -1;
    if (type == anaPot || (type == unknown && Worker->isAnaPotAv()))
    {
        int nPotCoeff, nLRCoeff, *pLRCoeff, n, r, NR, NAdCorr, TAdCorr = 0, PAdCorr = 0, NSpinRGamma;
        double Ri, Ra, iExp, iOff, iCoeff, Tm, De, *PotCoeff, *LRCoeff, A_ex = 0.0, beta = 0.0, gamma = 0.0, *adCorr, AdCorr_b, AdCorr_Rm;
        double RIso1AdCorr = 0.0, RIso2AdCorr = 0.0, *SpinRGamma, SpinRR1 = 0.0, SpinRR2 = 0.0, q_e = 0.0, q_f = 0.0;
        bool *LRCFree, *CoeffFree, ExcIntFree = false, *adCorrFree;
        QString Buffer;
        NR = Tab->rowCount();
        QTableWidgetItem *I;
        for (nPotCoeff = 0; (nPotCoeff < NR ? ((I = Tab->item(nPotCoeff, 0)) != 0 ? 
                I->text().left(1) == "a" : false) : false); nPotCoeff++) ;
        //printf("NR = %d, nPotCoeff=%d\n", NR, nPotCoeff);
        PotCoeff = new double[nPotCoeff];
        CoeffFree = new bool[nPotCoeff];
        for (n=0; n < nPotCoeff; n++) 
        {
            PotCoeff[n] = ((I = Tab->item(n, 1)) != 0 ? I->text().toDouble() : 0.0);
            CoeffFree[n] = (I != 0 ? I->icon().isNull() : true);
        }
        if (NR < n + 8) return;
        AdCorr_b = ((I = Tab->item(n, 1)) != 0 ? I->text().toDouble() : 0.0);
        AdCorr_Rm = ((I = Tab->item(n+1, 1)) != 0 ? I->text().toDouble() : 0.0);
        Tm = ((I = Tab->item(n+2, 1)) != 0 ? I->text().toDouble() : 0.0);
        Ri = ((I = Tab->item(n+3, 1)) != 0 ? I->text().toDouble() : 0.0);
        iExp = ((I = Tab->item(n+4, 1)) != 0 ? I->text().toDouble() : 0.0);
        iOff = ((I = Tab->item(n+5, 1)) != 0 ? I->text().toDouble() : 0.0);
        iCoeff = ((I = Tab->item(n+6, 1)) != 0 ? I->text().toDouble() : 0.0);
        Ra = ((I = Tab->item(n+7, 1)) != 0 ? I->text().toDouble() : 0.0);
        for (r=n+8; (r < NR ? ((I = Tab->item(r, 0)) != 0 ? I->text().left(1) == "C" : false) : false);
                 r++) ;
        LRCoeff = new double[nLRCoeff = r - n - 8];
        pLRCoeff = new int[nLRCoeff];
        LRCFree = new bool[nLRCoeff];
        LRCTabOffs = (r -= nLRCoeff);
        for (n=0; n < nLRCoeff; n++) 
        {
            Buffer = Tab->item(n+r, 0)->text();
            pLRCoeff[n] = Buffer.mid(1, Buffer.indexOf(" ") - 1).toInt();
            LRCoeff[n] = ((I = Tab->item(n+r, 1)) != 0 ? I->text().toDouble() : 0.0);
            LRCFree[n] = Tab->item(n+r, 0)->icon().isNull();
        }
        De = (NR > r + n ? ((I = Tab->item(n+r, 1)) != 0 ? I->text().toDouble() : 0.0) : 0.0);
        if (NR > (r+=n) + 3) if (Tab->item(r+1, 0)->text().left(4) == "A_ex")
        {
            ExcIntATabPos = (++r);
            A_ex = Tab->item(r, 1)->text().toDouble();
            ExcIntFree = Tab->item(r, 0)->icon().isNull();
            gamma = Tab->item(++r, 1)->text().toDouble();
            beta = Tab->item(++r, 1)->text().toDouble();
        }
        if (NR > r+1 ? Tab->item(r+1, 0)->text().left(3) == "q_e" : false) q_e = Tab->item(++r, 1)->text().toDouble();
        if (NR > r+1 ? Tab->item(r+1, 0)->text().left(3) == "q_f" : false) q_f = Tab->item(++r, 1)->text().toDouble();
        for (NAdCorr = 0, r++; (NR > r + NAdCorr ? Tab->item(r + NAdCorr, 0)->text().left(6) == "AdCorr" : false); NAdCorr++) ;
        if (NAdCorr > 0)
        {
            adCorr = new double[NAdCorr];
            adCorrFree = new bool[NAdCorr];
            adCorrTabOffs = n;
            for (n=0; n < NAdCorr; n++, r++) 
            {
                adCorrFree[n] = Tab->item(r, 0)->icon().isNull();
                adCorr[n] = (Tab->item(r, 1) != 0 ? Tab->item(r, 1)->text().toDouble() : 0.0);
            }
            TAdCorr = Tab->item(r++, 1)->text().toInt();
            PAdCorr = Tab->item(r++, 1)->text().toInt();
            RIso1AdCorr = Tab->item(r++, 1)->text().toDouble();
            RIso2AdCorr = Tab->item(r++, 1)->text().toDouble();
        }
        else 
        {
            adCorr = 0;
            adCorrFree = 0;
        }
        for (NSpinRGamma = 0; (NR > r + NSpinRGamma ? 
            Tab->item(r + NSpinRGamma, 0)->text().left(10) == "SpinRGamma" : false); NSpinRGamma++) ;
        if (NSpinRGamma > 0 && NR >= r + NSpinRGamma + 2)
        {
            
            SpinRGamma = new double[NSpinRGamma];
            for (n=0; n < NSpinRGamma; n++) SpinRGamma[n] = Tab->item(r+n, 1)->text().toDouble();
            SpinRR1 = Tab->item(r+n, 1)->text().toDouble();
            SpinRR2 = Tab->item(r+n+1, 1)->text().toDouble();
        }
        else SpinRGamma = 0;
        AnaPot *aPot;
        if (Worker->isAnaPotAv()) aPot = dynamic_cast<AnaPot*>(Worker);
        else setWorker(aPot = new AnaPot(Fit)); 
        aPot->setInnerWall(Ri, iOff, iExp, iCoeff);
        aPot->setPotCoeff(AdCorr_b, AdCorr_Rm, Tm, De, nPotCoeff, PotCoeff, CoeffFree);
        aPot->setLRCoeff(Ra, nLRCoeff, pLRCoeff, LRCoeff, LRCFree);
        if (A_ex != 0.0 || ExcIntFree) aPot->setExchangeCoeff(A_ex, beta, gamma, ExcIntFree);
        Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorrFree);
        Worker->setLambdaDoubling(q_e, q_f);
        Worker->setSpinRGamma(NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2, AdCorr_b, AdCorr_Rm);
    }
    else if (type == unknown ? (Tab->rowCount() > 6 ? (Tab->item(2, 0) != 0 ? 
                Tab->item(2, 0)->text().indexOf("R_ref") >= 0 : false) : false) : type == morseLongRange)
    {
        MLRPot *LPot = dynamic_cast<MLRPot*>(Worker);
        if (LPot == 0) setWorker(LPot = new MLRPot(Fit));
        QString B;
        int n, m, r, c, Nb, NCi, p = (Tab->item(0, 1) != 0 ? Tab->item(0, 1)->text().toInt() : 0);
        int q = (Tab->item(1, 1) != 0 ? Tab->item(1, 1)->text().toInt(): 0), NAdCorr;
        double Re = (Tab->item(3, 1) != 0 ? Tab->item(3, 1)->text().toDouble() : 0.0);
        double De = (Tab->item(4, 1) != 0 ? Tab->item(4, 1)->text().toDouble() : 0.0);
        double Rref = (Tab->item(2, 1) != 0 ? Tab->item(2, 1)->text().toDouble() : 0.0);
        double T = (Tab->item(5, 1) != 0 ? Tab->item(5, 1)->text().toDouble() : 0.0);
        for (n=6, Nb = 0; (n < Tab->rowCount() ? (Tab->item(n, 0) != 0 ? 
                Tab->item(n, 0)->text().indexOf("beta") != -1 : false) : false); n++, Nb++) ;
        for (NCi = 0; (n < Tab->rowCount() ? (Tab->item(n, 0) != 0 ?
                Tab->item(n, 0)->text().indexOf("C") != -1 : false) : false); n++, NCi++) ;
        double *b = new double[Nb], *Ci = new double[NCi];
        int *pCi = new int[NCi];
        for (n=0, r=6; n < Nb; n++, r++) 
            b[n] = (Tab->item(r, 1) != 0 ? Tab->item(r, 1)->text().toDouble() : 0.0);
        LRCTabOffs = r;
        for (n=0; n < NCi; n++, r++)
        {
            B = (Tab->item(r, 0) ? Tab->item(r, 0)->text() : "");
            c = B.indexOf('_');
            m = B.indexOf(' ', c+1);
            pCi[n] = B.mid(c+1, m - c - 1).toInt();
            Ci[n] = (Tab->item(r, 1) != 0 ? Tab->item(r, 1)->text().toDouble() : 0.0);
        }
        for (NAdCorr = 0, r++; (Tab->rowCount() > r + NAdCorr ? Tab->item(r + NAdCorr, 0)->text().left(6) == "AdCorr" : false); NAdCorr++) ;
        if (NAdCorr > 0)
        {
            double *adCorr = new double[NAdCorr];
            bool *adCorrFree = new bool[NAdCorr];
            adCorrTabOffs = r;
            for (n=0; n < NAdCorr; n++)
            {
                adCorrFree[n] = Tab->item(r, 0)->icon().isNull();
                adCorr[n] = Tab->item(r++, 1)->text().toDouble();
            }
            int TAdCorr = Tab->item(r++, 1)->text().toInt();
            int PAdCorr = Tab->item(r++, 1)->text().toInt();
            double RIso1AdCorr = Tab->item(r++, 1)->text().toDouble();
            double RIso2AdCorr = Tab->item(r++, 1)->text().toDouble();
            double AdCorr_b = Tab->item(r++, 1)->text().toDouble();
            double AdCorr_Rm = Tab->item(r++, 1)->text().toDouble();
            Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
        }
        else Worker->setAdCorr(0, 0, 0, 0, 0.0, 0.0);
        LPot->setCoefficients(p, q, T, De, Re, Rref, Nb, b, NCi, pCi, Ci);
    }
    else if (type == splinePot || Worker->isSplinePot() || (type == unknown && testIfSplinePot()))
    {
        double x, iA = 0.0, iO = 0.0, iExp = 6.0, Uinf = 0.0, RIso1AdCorr = 0.0, RIso2AdCorr = 0.0, AdCorr_b = 0.0, AdCorr_Rm = 0, Ro = 0.0, *LRC, *adCorr;
        int i, n, N = Tab->rowCount(), NAdCorr, numSplinePoints, NLRC, TAdCorr = 0, PAdCorr = 0, *PLRC;
        bool *LRCfree, *adCorrFree, calcYss = true, nanfound = false;
        QString Buffer[N];
        if (!Worker->isSplinePot()) setWorker(new SplinePot(Fit));
        Tab->setColumnCount(3);
        Tab->setHorizontalHeaderLabels(QStringList() << "R [A]" << "E [cm^-1]" << "d^2E/dR^2");
        for (n = NAdCorr = numSplinePoints = NLRC = 0; n<N; n++) 
        {
            Buffer[n] = Tab->item(n, 0)->text();
            if (Buffer[n].toDouble() > 0.0) numSplinePoints++;
            else if (Buffer[n].indexOf("C_") != -1) NLRC++;
            else if (Buffer[n].indexOf("iA") != -1) iA = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("iO") != -1) iO = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("iExp") != -1) iExp = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("U_inf") != -1) Uinf = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("Type AdCorr") != -1) TAdCorr = Tab->item(n, 1)->text().toInt();
            else if (Buffer[n].indexOf("Pow AdCorr") != -1) PAdCorr = Tab->item(n, 1)->text().toInt();
            else if (Buffer[n].indexOf("RIso1AdCorr") != -1) RIso1AdCorr = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("RIso2AdCorr") != -1) RIso2AdCorr = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("AdCorr b") != -1) AdCorr_b = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("AdCorr Rm") != -1) AdCorr_Rm = Tab->item(n, 1)->text().toDouble();
            else if (Buffer[n].indexOf("AdCorr") != -1) NAdCorr++;
            else if (Buffer[n].indexOf("R_o") != -1) Ro = Tab->item(n, 1)->text().toDouble();
        }
        SplinePoint *points = new SplinePoint[numSplinePoints];
        if (NLRC > 0)
        {
            LRC = new double[NLRC];
            PLRC = new int[NLRC];
            LRCfree = new bool[NLRC];
        }
        else
        {
            LRC = 0;
            PLRC = 0;
            LRCfree = 0;
        }
        if (NAdCorr > 0) 
        {
            adCorr = new double[NAdCorr];
            adCorrFree = new bool[NAdCorr];
        }
        else 
        {
            adCorr = 0;
            adCorrFree = 0;
        }
        for (n = numSplinePoints = NLRC = 0; n<N; n++)
        {
            if ((x = Buffer[n].toDouble()) > 0.0)
            {
                points[numSplinePoints].x = x;
                points[numSplinePoints].y = Tab->item(n, 1)->text().toDouble();
                if (points[numSplinePoints].y == 0.0)
                {
                    QLocale german(QLocale::German);
                    points[numSplinePoints].y = german.toDouble(Tab->item(n, 1)->text());
                }
                points[numSplinePoints].yss = (Tab->item(n, 2) != 0 ? Tab->item(n, 2)->text().toDouble() : 0.0);
                if (points[numSplinePoints].yss != 0.0) calcYss = false;
                if (isnan(points[numSplinePoints].yss) != 0 || isinf(points[numSplinePoints].yss) != 0) nanfound = true;
                points[numSplinePoints++].variable = Tab->item(n, 0)->icon().isNull();
            }
            else if ((i = Buffer[n].indexOf("C_")) != -1)
            {
                if (NLRC == 0) LRCTabOffs = n;
                PLRC[NLRC] = Buffer[n].right(Buffer[n].count() - i - 2).toInt();
                LRCfree[NLRC] = Tab->item(n, 0)->icon().isNull();
                LRC[NLRC++] = (Tab->item(n, 1) != 0 ? Tab->item(n, 1)->text().toDouble() : 0.0);
            }
            else for (i=0; i < NAdCorr; i++) if (Buffer[n] == "AdCorr" + QString::number(i)) 
            {
                if (i==0) adCorrTabOffs = n;
                adCorrFree[i] = Tab->item(n, 0)->icon().isNull();
                adCorr[i] = (Tab->item(n, 1) != 0 ? Tab->item(n, 1)->text().toDouble() : 0.0);
            }
        }
        if (nanfound) calcYss = true;
        if (NLRC == 0 && numSplinePoints > 2 ? 
                  points[numSplinePoints - 1].x - points[numSplinePoints - 2].x > 
                      10.0 * (points[numSplinePoints - 2].x - points[numSplinePoints - 3].x) : false)
        {
            Tab->blockSignals(true);
            Uinf = points[numSplinePoints - 1].y;
            Tab->item(--numSplinePoints, 0)->setText("U_inf");
            Tab->item(numSplinePoints, 2)->setText("");
            Tab->blockSignals(false);
        }
        Worker->setSplinePot(numSplinePoints, points, NLRC, PLRC, LRC, iA, iO, iExp);
        Worker->setLRCFree(LRCfree);
        Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorrFree);
        Worker->setUinf(Uinf);
        Worker->setR_o(Ro);
        if (calcYss) calcyss(false);
    }
    else if (type == unknown ? Tab->item(0, 0)->text().indexOf("A0") >= 0 || Tab->item(0, 0)->text().indexOf("B") >= 0 : type == modTangToen)
    {
        double *A, *C, alpha, beta, gamma, Uinf, *SA = 0, *SC = 0, Salpha = 0.0;
        double Sinf = 0.0, Sbeta = 0.0, B = 0.0, alpha1 = 0.0, *adCorr, RIso1AdCorr = 0.0, RIso2AdCorr = 0.0, AdCorr_Rm = 0.0, AdCorr_b = 0.0;
        int NA, NC, n=0, m, *pC, NSA = 0, NSC = 0, *pSC = 0, NAdCorr, TAdCorr = 0, PAdCorr = 0;
        bool *adCorrFree;
        QString Buffer;
        if (Tab->item(0, 0)->text().indexOf("B") >= 0) B = Tab->item(n++, 1)->text().toDouble();
        if (Tab->item(n, 0)->text().indexOf("alpha1", 0, Qt::CaseInsensitive) >= 0) 
            alpha1 = Tab->item(n++, 1)->text().toDouble();
        for (NA = 0; ((m = Tab->item(n + NA, 0)->text().indexOf("A")) >= 0 ? 
            Tab->item(n + NA, 0)->text().mid(m+1, 1).toInt() == NA : false); NA++) ;
        if (Tab->item(n + NA, 0)->text().indexOf("alpha", 0, Qt::CaseInsensitive) >= 0
            && Tab->item(n + NA + 1, 0)->text().indexOf("beta", 0, Qt::CaseInsensitive) >= 0
            /*&& Tab->item(NA + 2, 0)->text().indexOf("gamma", 0, Qt::CaseInsensitive) >= 0
            && Tab->item(NA + 3, 0)->text().indexOf("C") >= 0*/)
        {
            for (NC = 1; Tab->item(n + NA + 3 + NC, 0)->text().indexOf("C") >= 0; NC++) ;
            if (Tab->item(n + NA + 3 + NC, 0)->text().indexOf("U_inf") >= 0
                && !Tab->item(n + NA + 3 + NC, 1)->text().isEmpty())
            {
                A = new double[NA];
                pC = new int[NC];
                C = new double[NC];
                for (m=0; m < NA; m++) A[m] = Tab->item(n++, 1)->text().toDouble();
                alpha = Tab->item(n++, 1)->text().toDouble();
                beta = Tab->item(n++, 1)->text().toDouble();
                if (Tab->item(n, 0)->text().indexOf("gamma", 0, Qt::CaseInsensitive) >= 0)
                    gamma = Tab->item(n++, 1)->text().toDouble();
                else gamma = 0;
                LRCTabOffs = n;
                for (m=0; m < NC; m++)
                {
                    Buffer = Tab->item(n, 0)->text();
                    pC[m] = Buffer.right(Buffer.length() 
                                            - Buffer.indexOf("C") - 1).toInt();
                    C[m] = Tab->item(n++, 1)->text().toDouble();
                }
                Uinf = Tab->item(n++, 1)->text().toDouble();
                if (Tab->rowCount() >= n+5)
                {
                    while ((m = Tab->item(n + NSA, 0)->text().indexOf("SA")) >= 0 ? 
                           Tab->item(n + NSA, 0)->text().mid(m+2, 1).toInt() == NSA : false)
                        NSA++;
                    if (NSA > 0)
                    {
                        SA = new double[NSA];
                        for (m=0; m < NSA; m++) SA[m] = Tab->item(n++, 1)->text().toDouble();
                    }
                    if (Tab->item(n, 0)->text().indexOf("Salpha", 0, Qt::CaseInsensitive) >= 0)
                        Salpha = Tab->item(n++, 1)->text().toDouble();
                    if (Tab->item(n, 0)->text().indexOf("Sbeta", 0, Qt::CaseInsensitive) >= 0)
                        Sbeta = Tab->item(n++, 1)->text().toDouble();
                    while (Tab->item(n + NSC, 0)->text().indexOf("SC", 0, Qt::CaseInsensitive) >= 0)
                        NSC++;
                    if (NSC > 0)
                    {
                        pSC = new int[NSC];
                        SC = new double[NSC];
                        for(m=0; m < NSC; m++)
                        {
                            Buffer = Tab->item(n, 0)->text();
                            pSC[m] = Buffer.right(Buffer.length() 
                                            - Buffer.indexOf("C") - 1).toInt();
                            SC[m] = Tab->item(n++, 1)->text().toDouble();
                        }
                    }
                    if (Tab->item(n, 0)->text().indexOf("S_inf", 0, Qt::CaseInsensitive) >= 0)
                        Sinf = Tab->item(n++, 1)->text().toDouble();
                }
                for (NAdCorr = 0; (Tab->rowCount() > n + NAdCorr + 6 ? Tab->item(n + NAdCorr, 0)->text().left(6) == "AdCorr" : false); NAdCorr++) ;
                if (NAdCorr > 0)
                {
                    adCorr = new double[NAdCorr];
                    adCorrFree = new bool[NAdCorr];
                    adCorrTabOffs = n;
                    for (m=0; m < NAdCorr; m++) 
                    {
                        adCorrFree[m] = Tab->item(n, 0)->icon().isNull();
                        adCorr[m] = Tab->item(n++, 1)->text().toDouble();
                    }
                    TAdCorr = Tab->item(n++, 1)->text().toInt();
                    PAdCorr = Tab->item(n++, 1)->text().toInt();
                    RIso1AdCorr = Tab->item(n++, 1)->text().toDouble();
                    RIso2AdCorr = Tab->item(n++, 1)->text().toDouble();
                    AdCorr_b = Tab->item(n++, 1)->text().toDouble();
                    AdCorr_Rm = Tab->item(n++, 1)->text().toDouble();
                }
                else 
                {
                    adCorr = 0;
                    adCorrFree = 0;
                }
                MTTPot *mTTPot = new MTTPot(Fit);
                mTTPot->setCoefficients(B, alpha1, NA, A, alpha, beta, gamma, NC, pC, C, Uinf, 
                                        NSA, SA, Salpha, Sbeta, NSC, pSC, SC, Sinf);
                setWorker(mTTPot);
                Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorrFree);
            }
        }
    }
    else if (type == unknown ? Tab->item(0, 0)->text()[0] == 'A' && Tab->rowCount() > 5 : type == tangToennies)
    {
        int n, m, N, NA, NR = Tab->rowCount(), NAdCorr, TAdCorr = 0, PAdCorr = 0;
        double *A, b, Uinf, *C, *adCorr, RIso1AdCorr = 0.0, RIso2AdCorr = 0.0, AdCorr_b = 0.0, AdCorr_Rm = 0.0;
        bool *adCorrFree;
        for (NA = 0; (NA < NR ? Tab->item(NA, 0)->text().left(1) == "A" : false); NA++) ;
        for (n=0, A = new double[NA]; n < NA; n++) A[n] = Tab->item(n, 1)->text().toDouble();
        b = Tab->item(NA, 1)->text().toDouble();
        for (N=0; (N + NA + 1 < NR ? Tab->item(N + NA + 1, 0)->text().left(2) == "C_" : false) ; N++) ;
        LRCTabOffs = NA + 1;
        for (n=0, C = new double[N]; n<N; n++) C[n] = Tab->item(n + LRCTabOffs, 1)->text().toDouble();
        if (N + NA + 1 < NR) Uinf = Tab->item(n + NA + 1, 1)->text().toDouble();
        else Uinf = 0.0;
        for (NAdCorr = 0; (NR > n + NA + NAdCorr + 8 ? 
            Tab->item(n + NA + NAdCorr + 2, 0)->text().left(6) == "AdCorr" : false); NAdCorr++) ;
        if (NAdCorr > 0)
        {
            adCorr = new double[NAdCorr];
            adCorrFree = new bool[NAdCorr];
            adCorrTabOffs = n;
            for (m=0, n += (NA + 2); m < NAdCorr; m++) 
            {
                adCorrFree[m] = Tab->item(n, 0)->icon().isNull();
                adCorr[m] = Tab->item(n++, 1)->text().toDouble();
            }
            TAdCorr = Tab->item(n++, 1)->text().toInt();
            PAdCorr = Tab->item(n++, 1)->text().toInt();
            RIso1AdCorr = Tab->item(n++, 1)->text().toDouble();
            RIso2AdCorr = Tab->item(n++, 1)->text().toDouble();
            AdCorr_b = Tab->item(n++, 1)->text().toDouble();
            AdCorr_Rm = Tab->item(n++, 1)->text().toDouble();
        }
        else 
        {
            adCorr = 0;
            adCorrFree = 0;
        }
        TangToenniesPot *TTPot = new TangToenniesPot(Fit);
        TTPot->setPotCoeff(A, NA, b, C, N, Uinf);
        setWorker(TTPot);
        Worker->setAdCorr(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorrFree);
    }
    else printf("Potential::UpdatePot error: unknown type=%d\n", type);
    Changed();
}

double Potential::getMinR()
{
    //if (aPot == 0 && LPot == 0 && mTTPot == 0 && points != 0) return points[0].x;
    return rmin;
}

double Potential::getMaxR()
{
    //if (aPot == 0 && LPot == 0 && points != 0) return points[numSplinePoints - 1].x;
    return rmax;
}

int Potential::getNFC()
{
    if (Worker->isFSA()) 
        return (State != 0 ? int(2.0 * State->getS()) + 1 : 1);
    return 1;
}

double* Potential::getAdCorrEnergy(double Rmin, double Step, int numPoints, double mIso1, double mIso2)
{
    double *R = new double[numPoints];
    int n;
    for (n=0; n < numPoints; n++) R[n] = 0.0;
    Worker->CalcAdCorrPot(R, R, Rmin, Step, numPoints, mIso1, mIso2);
    return R;
}

void Potential::getReDe(double &Re, double &De)
{
    double Te;
    getMinimum(Re, Te);
    if (Worker->isLPot()) De = dynamic_cast<MLRPot*>(Worker)->GetDe();
    else if (Worker->isTTPot()) De = -Te;
    else 
    {
        double Uinf = Worker->getUinf();
        De = Uinf - Te;
    }
}

double Potential::getMinimumE()
{
    double R, E;
    getMinimum(R, E);
    return E;
}

double Potential::guessBe()
{
    if (molecule == 0) return -1;
    IsoTab *Iso = molecule->getIso();
    double Re, E;
    getMinimum(Re, E);
    double Be = 1e18 * C_h / (M_PI * M_PI * 8.0 * Iso->redMass[Iso->refIso] * Re * Re * C_u * C_c);
    delete Iso;
    return Be;
}

void Potential::setMinimum(double E)
{
    shiftEnergyOffset(E - getMinimumE());
}

double Potential::getAsymptote()
{
    return Worker->getUinf();
    /*if (points != 0)
    {
        int n = numSplinePoints - 1;
        if (SLIC6) return points[n].y - points[n].yss * pow(points[n].x, 6.0);
        return points[n].y;
        return Uinf;
    }
    printf("Potential::getasymptote: error, there is no potential data available!");
    return 0.0;*/
}

void Potential::setCoefficients()
{
    Tab->blockSignals(true);
    QStringList Set;
    int n, N = Tab->rowCount();
    mWasMoving = false;
    for (n=0; n<N; n++) Set << Tab->item(n, 0)->text();
    CoefficientDialog *D = new CoefficientDialog(this, Set);
    if (D->exec() == QDialog::Rejected) return;
    Set = D->getResults();
    if (Set.count() > 6 ? Set[2].indexOf("R_ref") == -1 : true) Set.sort();
    Tab->setRowCount(N = Set.size());
    for (n=0; n<N; n++) Tab->setItem(n, 0, new QTableWidgetItem(Set[n]));
    Tab->blockSignals(false);
    updatePot();
    delete D;
}

void Potential::setAsymptote(double E)
{
    shiftEnergyOffset(E - getAsymptote());
}

void Potential::shiftEnergyOffset(double E)
{
    Worker->shift(E);
    UpdateTab();
    Changed();
}

ElState *Potential::getElState()
{
    return State;
}

void Potential::autoCalcScatLengthsPotentialSet(int NumWFPoints)
{
    IsoTab *IsoT = Worker->getIsoT();
    int i, c, n, m, nC, *lExp, NIso = IsoT->numIso;
    double Ri, Ra, *lC, Rmin, Emin, iO, iC, Rm, Tm, b, iExp;
    AnaPot *aPot;
    QString PotFile, IsoS1[NIso], IsoS2[NIso], Buffer, IBuff;
    QString PotDir = QFileDialog::getExistingDirectory(this, 
            "Please select the directory of potentials to use");
    if (PotDir.isEmpty()) return;
    PotDir.replace('/', '\\');
    QString ProgDir = QFileDialog::getOpenFileName(this,
            "Please give the path of the program 'scat_all_asym9.exe'", PotDir, "Programs (*.exe)");
    ProgDir.replace('/', '\\');
    if (ProgDir.isEmpty()) return;
    QString InFile = QFileDialog::getOpenFileName(this,
            "Please give the input file for the program", ProgDir);
    if (InFile.isEmpty()) return;
    InFile.replace('/', '\\');
    QString OutFile = InFile + ".out";
    QString ResFile = QFileDialog::getSaveFileName(this, 
            "Please select the file name for the results", InFile);
    if (ResFile.isEmpty()) return;
    ResFile.replace('/', '\\');
    QFile IFile(InFile);
    QFile OFile(OutFile);
    IFile.open(QIODevice::ReadOnly);
    QTextStream IS(&IFile);
    QTextStream OS(&OFile);
    QStringList InList;
    while ((IBuff = IS.readLine()).left(13) != "  AtomMass_A=") InList << IBuff;
    for (n=0; n<3; n++) IS.readLine();
    IBuff = IS.readAll();
    IFile.close();
    QFile RFile(ResFile);
    QProcess *Calc = new QProcess;
    Calc->setWorkingDirectory(InFile.left(InFile.lastIndexOf('\\')));
    RFile.open(QIODevice::WriteOnly);
    QTextStream RS(&RFile);
    QDir PDir(PotDir);
    QFileInfoList PList = PDir.entryInfoList(QDir::Dirs | QDir::Files | QDir::NoDotAndDotDot);
    RS << "Potential\tC6\tC8\tC10\tRmin\tDe\tD0\tn\tb\tRm\tTm\tnC\tRi\tRa";
    for (n=0; n < NIso; n++) 
    {
        RS << "\t" << IsoT->getIsoName(n);
        IsoS1[n] = "  AtomMass_A=" + QString::number(IsoT->mIso1[n], 'f', 8) + "\n";
        IsoS2[n] = "  AtomMass_B=" + QString::number(IsoT->mIso2[n], 'f', 8) + "\n";
    }
    RS << "\n";
    for (n=0; n < PList.size(); n++) 
    {
        if (PList[n].isDir())
        {
            PDir.setPath(PList[n].filePath());
            PList << PDir.entryInfoList(QDir::Dirs | QDir::Files | QDir::NoDotAndDotDot);
            continue;
        }
        if (!PList[n].isFile()) continue;
        if (!readData(PotFile = PList[n].filePath())) continue;
        aPot = Worker->getAnaPotForReading();
        if (aPot == 0) continue;
        if (!isImported()) writePotData(PotFile = ProgDir.left(ProgDir.lastIndexOf('\\') + 1) + "TempPot.dat");
        aPot->getLRCoeffForReading(nC, lExp, lC, &Ra);
        RS << PList[n].fileName() << "\t" << (nC > 0 ? QString::number(lC[0], 'f', 8) : "0.0") << "\t";
        RS << (nC > 1 ? QString::number(lC[1], 'f', 8) : "0.0");
        RS << "\t" << (nC > 2 ? QString::number(lC[2], 'f', 8) : "0.0") << "\t";
        Worker->unlock();
        aPot->getMinimum(Rmin, Emin);
        RS << QString::number(Rmin, 'f', 8) << "\t" << QString::number(aPot->getUinf() - Emin, 'f', 5);
        RS << "\t" << QString::number(getD0(NumWFPoints), 'f', 5) << "\t";
        aPot->getInnerWall(Ri, iO, iExp, iC);
        RS << QString::number(iExp, 'g', 10) << "\t";
        aPot->getPotCoeff(b, Rm, Tm, Emin, nC, lC);
        RS << QString::number(b, 'f', 2) << "\t" << QString::number(Rm, 'f', 8) << "\t";
        RS << QString::number(Tm, 'f', 6) << "\t" << nC << "\t";
        RS << QString::number(Ri, 'f', 3) << "\t" << QString::number(Ra, 'f', 6);
        for (m=0; m < NIso; m++)
        {
            IFile.open(QIODevice::WriteOnly);
            IS.setDevice(&IFile);
            for (i=0; i < InList.size(); i++) IS << InList[i] << "\n";
            IS << IsoS1[m] << IsoS2[m] << " /\n" << "   11" << PotFile.replace('/', '\\') << "\n" 
                    << IBuff;
            IFile.close();
            Calc->start(ProgDir);
            if (!Calc->waitForStarted(5000))
            {
                QMessageBox::information(this, "MolSpektAnalysis", "0Error starting " + ProgDir);
                return;
            }
            for (c=0, Buffer = ""; c < 10 && Buffer.isEmpty(); Buffer = Calc->readLine(), c++)
            {
                if (!Calc->waitForReadyRead(5000))
                {
                    QMessageBox::information(this, "MolSpektAnalysis", "1Error starting " + ProgDir);
                    return;
                }
            }
            if (!Buffer.contains("input:"))
            {
                QMessageBox::information(this, "MolSpektAnlalyisis", 
                                         "2Error in " +  ProgDir + ": " + Buffer);
                return;
            }
            Calc->write((InFile + "\n").toUtf8());
            if (!Calc->waitForReadyRead(5000))
            {
                QMessageBox::information(this, "MolSpektAnalysis", "3Error starting " + ProgDir);
                return;
            }
            Buffer = Calc->readAll();
            if (!Buffer.contains("output:"))
            {
                QMessageBox::information(this, "MolSpektAnlalyisis", 
                                         "4Error in " +  ProgDir + ": " + Buffer);
                return;
            }
            Calc->write((OutFile + "\n").toUtf8());
            /*Calc->waitForReadyRead(5000);
            QMessageBox::information(this, "MolSpektAnalysis", 
                                         "Result in " + ProgDir + ": " + Calc->readAll());
            Calc->waitForReadyRead(5000);
            Calc->readAll();*/
            if (!Calc->waitForFinished(5000))
            {
                QMessageBox::information(this, "MolSpektAnalysis", 
                                         "5Error in " + ProgDir + ": " + Calc->readAll());
                return;
            }
            OFile.open(QIODevice::ReadOnly);
            OS.setDevice(&OFile);
            while ((Buffer = OS.readLine()).left(48) != 
                             " Elastic and Inelastic Scattering Lengths [bohr]" && !OS.atEnd()) ;
            if (Buffer.indexOf("Scattering Lengths") == -1) 
            {
                QMessageBox::information(this, "MolSpektAnalysis", "1Error reading output file, InFile=" + InFile + ", OutFile=" + OutFile);
                return;
            }
            while ((Buffer = OS.readLine()).left(18) != "  | 0  0>  | 0  0>") ;
            if (Buffer.left(18) != "  | 0  0>  | 0  0>")
            {
                QMessageBox::information(this, "MolSpektAnalysis", "2Error reading output file");
                return;
            }
            OFile.close();
            RS << "\t" << Buffer.right(Buffer.length() - 18);
        }
        RS << "\n";
    }
    delete Calc;
}

void Potential::addPoint(double x, double y)
{
    if (Worker->addPoint(x, y)) calcyss(mWasMoving = true);
    else QMessageBox::information(this, "MolSpektAnalysis", "The spline points cannot be changed manually while a fit is running.");
}

void Potential::movePoint(int &n, double x, double y)
{
    if (Worker->movePoint(n, x, y))
    {
        Worker->calcIWallByTwoSplinePoints();
        calcyss(mWasMoving = true);
    }
    else QMessageBox::information(this, "MolSpektAnalysis", "The spline points cannot be changed manually while a fit is running.");
}

void Potential::removePoint(int n)
{
    if (Worker->removePoint(n)) calcyss(mWasMoving = true);
    else QMessageBox::information(this, "MolSpektAnalysis", "The spline points cannot be changed manually while a fit is running.");
}

void Potential::createMLRPot()
{
    //if (LPot != 0) {LPot->test(); return;}
/*    int n, m, nLRC, *pLRC;
    double *lRC;
    QStringList L;
    int p, q, Nbeta;
    double TAs, De, Re, Rref, *beta;
    getLRCoeffForReading(nLRC, pLRC, lRC);
    if (Worker->isLPot()) dynamic_cast<MLRPot*>(Worker)->getCoefficients(p, q, TAs, De, Re, Rref, Nbeta, beta, nLRC, pLRC, lRC);
    for (n=0; n < nLRC; n++) 
    {
        switch (pLRC[n])
        {
            case 6:
                if (L.size() == 0) L << "";
                break;
            case 8:
                for (m = L.size(); m<2; m++) L << "";
                break;
            case 10:
                for (m = L.size(); m<3; m++) L << "";
                break;
            case 12:
                for (m = L.size(); m<4; m++) L << "";
                break;
        }
        L << QString::number(lRC[n], 'e', 6);
    }
    CMLRPotDialog *D = new CMLRPotDialog(&L, this);
    for (n=0; n < nLRC; n++) D->LRC->addItem("C" + QString::number(pLRC[n]) + " = " + QString::number(lRC[n], 'e', 6));
    if (Worker->isLPot())
    {
        D->pE->setText(QString::number(p));
        D->qE->setText(QString::number(q));
        D->RrefE->setText(QString::number(Rref, 'f', 1));
        D->NCE->setText(QString::number(Nbeta));
        delete[] beta;
    }
    if (pLRC != 0)
    {
        delete[] pLRC;
        delete[] lRC;
    }
    if (D->exec() == QDialog::Rejected) 
    {
        delete D;
        return;
    }
    double *R, mR = D->mRE->text().toDouble(), MR = D->MRE->text().toDouble(), r;
    int NP = D->NPE->text().toInt(), NC = D->NCE->text().toInt(), i, o, *PLRC;
    double h = (MR - mR) / double(NP - 1), EMin, UInf;
    nLRC = D->LRC->count();
//    MrqLev ML(NC + nLRC + 3, NP);
//    if ((ML.a[0] = UInf = getAsymptote()) != 0.0) ML.ia[0] = false;
    //ML.ia[0] = ML.ia[1] = ML.ia[2] = false;
//    getMinimum(ML.a[2], EMin);
//    getReDe(ML.a[2], ML.a[1]);
    QString S;
    R = getPoints(mR, MR, NP);
//    fitBeta(R, mR, MR, NP, EMin, beta = new double[NC], NC, p = D->pE->text().toInt(), q = D->qE->text().toInt(), ML.a[1], ML.a[2],
//            Rref = D->RrefE->text().toDouble(), lRC, pLRC, nLRC);
    for (n=0, m=3; n < NC; n++, m++) ML.a[m] = beta[n];
    for (pLRC = new int[nLRC], PLRC = new int[nLRC], lRC = new double[nLRC], n=0; n < nLRC; n++) 
    {
        i = ((i = (S = D->LRC->item(n)->text()).indexOf(' ')) > 0 ? i : S.length());
        PLRC[n] = pLRC[n] = S.mid(1, i-1).toInt();
        ML.a[m] = lRC[n] = ((i = S.indexOf('=')) > 0 ? S.right(S.length() - i - 2).toDouble() : 0.0);
        ML.ia[m++] = (lRC[n] == 0.0 ? true : false);
    }
    Potential *Pot = MW->CreatePotential();
    MLRPot *MP = new MLRPot(Pot->Fit);
    MP->setCoefficients(p, q, ML.a[0], ML.a[1], ML.a[2], Rref, NC, beta, nLRC, pLRC, lRC);
    //for (r = 4.0; r <= 10.0; r+=100*h) printf("R=%f, E=%e\n", r, MP->GetPoint(r)); return;
    beta = new double[NC];
    lRC = new double[nLRC];
    for (n=0, r = mR; n < NP; n++, r+=h)
    {
        ML.y[n] = R[n];
        ML.sig[n] = 0.01;
        ML.dyda[n][0] = 1.0;
        MP->getDerivatives(r, ML.E[n], beta, lRC, ML.dyda[n][2], ML.dyda[n][1]);
        for (m=0, o=3; m < NC; m++) ML.dyda[n][o++] = beta[m];
        for (m=0; m < nLRC; m++) ML.dyda[n][o++] = lRC[m];
        if (isinf(ML.E[n])) 
        {
            ML.E[n] = 0.0;
            for (m=0; m < NC + nLRC + 3; m++) ML.dyda[n][m] = 0.0;
        }
    }
    ML.init();
    printf("init: %e ma=%d, mfit=%d\n", ML.chisq, ML.ma, ML.mfit);
    //return;
    //fitMLRPot(mR, MR, R, NP, MP, fLRC);
    for (i=0; i < 1000; i++)
    {
        //for (n=0; n < ML.ma; n++) printf("%e ", ML.atry[n]); printf("\n");
        pLRC = new int[nLRC];
        for (m=0, o=3; m < NC; m++) beta[m] = ML.atry[o++];
        for (m=0; m < nLRC; m++)
        {
            pLRC[m] = PLRC[m];
            lRC[m] = ML.atry[o++];
        }
        MP->setCoefficients(p, q, ML.atry[0], ML.atry[1], ML.atry[2], Rref, NC, beta, nLRC, pLRC, lRC);
        beta = new double[NC];
        lRC = new double[nLRC];
        for (n=0, r = mR; n < NP; n++, r+=h)
        {
            MP->getDerivatives(r, ML.E[n], beta, lRC, ML.dyda[n][2], ML.dyda[n][1]);
            for (m=0, o=3; m < NC; m++) ML.dyda[n][o++] = beta[m];
            for (m=0; m < nLRC; m++) ML.dyda[n][o++] = lRC[m];
        }
        if (!ML.fit()) 
        {
            printf("finished: chisq=%e\n", ML.chisq);
            break;
        }
        if (ML.alamda > 1e300)
        {
            printf("lamda to large!\n");
            break;
        }
        printf("%d: %e nd=%d, l=%e\n", i, ML.chisq, ML.done, ML.alamda);
    }
    pLRC = new int[nLRC];
    for (m=0, o=3; m < NC; m++) beta[m] = ML.a[o++];
    for (m=0; m < nLRC; m++)
    {
        pLRC[m] = PLRC[m];
        lRC[m] = ML.a[o++];
    }
    MP->setCoefficients(p, q, UInf, ML.a[1], ML.a[2], Rref, NC, beta, nLRC, pLRC, lRC);
    Pot->setWorker(MP);
    Pot->UpdateTab();
    State->addPotential(Pot);
    Pot->show();
    delete[] R;
    delete[] PLRC;
    delete D;*/
}

PotentialType Potential::getPotType()
{
    if (Worker->isSplinePot()) return SplinePotential;
    if (Worker->isAnaPotAv()) return analyticalPotential;
    if (Worker->isTTPot()) return TangToenniesPotential;
    if (Worker->isLPot()) return MorseLongRangePotential;
    if (Worker->isMTTPot()) return ModifiedTangToenniesPotential;
    return NoPotential;
}

FitResult Potential::getFitResult()
{
    FitResult Result = Worker->getFitResult();
    Result.ProcessNum = threadNum;
    return Result;
}
