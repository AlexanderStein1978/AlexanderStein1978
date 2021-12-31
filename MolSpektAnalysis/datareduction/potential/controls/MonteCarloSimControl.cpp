//
// C++ Implementation: MonteCarloSimControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "MonteCarloSimControl.h"
#include "potential.h"
#include "elstate.h"
#include "fitdata.h"
#include "utils.h"
#include "MainWindow.h"


MonteCarloSimControl::MonteCarloSimControl(MainWindow* mw, Potential* Pot, QString PotDir, int NumIterations, int NumParFits, 
                                           double UncFact, int NumWFPoints) : ThreadControl(mw, Pot->getElState(), NumIterations, NumParFits)
{
    int n;
    double X, S, SF = 32768 / 2506.62827;
    QStringList HL;
    setWindowTitle("MCS Control");
    UncF = UncFact;
    Pot->getFQS(NumWFPoints);
    MFD = St->getFitData();
    MCSFD = MW->CreateFitData();
    MCSFD->setName("MCSFD_" + Pot->getName());
    MCSFD->setSource(MFD->getName() + " from " + MFD->getSource());
    MCSFD->setvMax(MFD->getvMax());
    MCSFD->setJMax(MFD->getJMax());
    St->setMainFitData(MCSFD);
    PotName = Pot->getName();
    PotFileName = Pot->getFName();
    ResultDir = PotDir;
    NumPDRows = NumPDColumns = 0;
    ResultTab->setWindowTitle("Results Monte Carlo sim " + Pot->getName());
    PotData = Pot->getData(NumPDRows, NumPDColumns);
    Pot->getFixedCoefficients(NFC, FCRows);
    ResultTab->setTabDimensions(NumIt, NumResultTabColumns = NumPDRows + 1);
    for (n=0; n < NumPDRows; n++) HL << PotData[n][0];
    HL << "FQS";
    ResultTab->setHorizontalHeader(HL);
    MCSTab = new double[32768];
    for (n=1, X = -10.0, S = 0.0; X<=0.0; X += 0.001)
    {
        S += SF * exp(-0.5*X*X);
        while (int(S) < n)
        {
            X += 0.001;
            S += SF * exp(-0.5*X*X);
        }
        while (n <= int(S)) MCSTab[n++] = X;
    }
    for (MCSTab[0] = MCSTab[1] - 1.0 ; n < 32768; n++) MCSTab[n] = -MCSTab[32767 - n];
    StartCalcs();
}

MonteCarloSimControl::~MonteCarloSimControl()
{
    Destroy(PotData, NumPDRows);
    delete[] MCSTab;
}

void MonteCarloSimControl::FitFinished(int ThreadNum, double FQS)
{
    QString File = ResultDir + PotFileName.left(PotFileName.lastIndexOf('.')) + "MCS" + QString::number(PotI[ThreadNum]) + ".pot";
    int NumR = 0, NumC = 0, n;
    QString **ResData = Pots[ThreadNum]->getData(NumR, NumC), RData[NumPDRows + 1];
    for (n=0; n < NumR; n++) RData[n] = ResData[n][1];
    RData[n] = QString::number(FQS, 'f', 6);
    ResultTab->setRowData(PotI[ThreadNum], RData);
    Pots[ThreadNum]->writeData(File);
    Destroy(ResData, NumR);
    CalcFinished(ThreadNum);
}

void MonteCarloSimControl::StartCalc(int p, int n)
{
    int NL; // i=0, j=0;
    TableLine *Lines;
    if (Pots[p] == 0)
    {
        Pots[p] = MW->CreatePotential(p);
        if (Pots[p] == 0) return;
        St->addPotential(Pots[p]);
        Pots[p]->show();
        connect(Pots[p], SIGNAL(FitFinished(int,double)), this, SLOT(FitFinished(int,double)));
        connect(Pots[p], SIGNAL(FitTerminated(int)), this, SLOT(CalcTerminated(int)));
    }
    Pots[p]->setName(PotName + "MCS" + QString::number(n));
        /*FileName = File + QString::number(n) + ".pot";
        Test.setFileName(FileName);
        if (Test.exists())
        {
            Pots[p]->readData(FileName);
            if (FD == 0) 
            {
                FD = new FitData;
                St->setMainFitData(FD);
            }
            FD->readData(File + QString::number(n) + ".fdat");
            Pots[p]->CalcFQS();
            FD->getData(Lines, NL);
            for (i=0, FQS = 0.0; i < NL; i++) FQS += Lines[i].DevR * Lines[i].DevR;
            delete[] Lines;
            FQS = Pots[p]->getFQS();
            QString RData[NumPDRows + 1], **ResData = Pots[p]->getData(i, j);
            for (i=0; i < NumPDRows; i++) RData[i] = ResData[i][1];
            RData[i] = QString::number(FQS, 'f', 6);
            Destroy(ResData, i);
            ResultTab->setRowData(n, RData);
            Progress->setValue(++ItFinished);
            continue;
        }
        else 
        {
            St->setMainFitData(MFDN);*/
    MFD->getData(Lines, NL);
    MCSFD->setData(Lines, NL);
    Pots[p]->setData(PotData, NumPDRows, NumPDColumns);
    if (NFC > 0) Pots[p]->VaryCoefficients(false, NFC, FCRows);
    Pots[p]->MonteCarloSimIteration(MCSTab, UncF);
}
