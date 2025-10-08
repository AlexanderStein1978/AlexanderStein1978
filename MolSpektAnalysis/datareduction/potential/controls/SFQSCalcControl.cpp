//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "SFQSCalcControl.h"
#include "potential.h"
#include "fitdata.h"
#include "utils.h"
#include "isotab.h"
#include "elstate.h"
#include "MainWindow.h"

#include <QDir>


SFQSCalcControl::SFQSCalcControl(MainWindow* MW, ElState* St, Potential *Pot, QString PotDir, QString FDatDir, double SFQSRad, int NumParFits)
                    : ThreadControl(MW, St, 0, NumParFits)
{
    QDir Dir(PDir = PotDir);
    int NRows = 0, NCols = 0, c;
    QString **PotData = Pot->getData(NRows, NCols);
    QStringList HeaderLabels;
    setWindowTitle("FQS Calc Control");
    OPot = Pot;
    FDir = FDatDir;
    Rad = SFQSRad;
    SFQSU = 0;
    AdTabRows = 1;
    if (PDir.right(1) != "\\" && PDir.right(1) != "/") PDir += DIRSEP;
    PotList = Dir.entryList(QStringList("*.pot"), QDir::Files);
    if (!FDir.isEmpty())
    {
        if (FDir.right(1) != "\\" && FDir.right(1) != "/") FDir += DIRSEP;
        MFD = MW->CreateFitData();
        MFD->setName("MCSFDat");
    }
    ResultTab->setWindowTitle("Results Monte Carlo sim " + Pot->getName());
    for (c=0; c < NRows; c++) HeaderLabels << PotData[c][0];
    HeaderLabels << "SFQS" << "FQS";
    ResultTab->setTabDimensions(PotList.count() + 1, NumResultTabColumns = NRows + 2);
    ResultTab->setHorizontalHeader(HeaderLabels);
    Destroy(PotData, NRows);
    connect(Pot, SIGNAL(SFQScalculated(int,double,double)), this, SLOT(CalcFinished(int,double,double)));
    Pot->calcFQS_SFQS(0, SFQSRad);
}

SFQSCalcControl::~SFQSCalcControl()
{
    int I, c, J, NI = OPot->getIsoT()->numIso, NC = int(St->getS() * 2.0) + 1;
    if (SFQSU != 0)
    {
        for (c=0; c < NC; c++)
        {
            for (I=0; I < NI; I++)
            {
                for (J=0; J < mJ[I]; J++) delete[] SFQSU[c][I][J];
                delete[] SFQSU[c][I];
            }
            delete[] SFQSU[c];
        }
        delete[] SFQSU;
        for (I=0; I < NI; I++) delete[] mv[I];
        delete[] mv;
        delete[] mJ;
    }
}

void SFQSCalcControl::CalcFinished(int ThreadNum, double SFQS, double FQS)
{
    int n, r = (ThreadNum >= 0 ? PotI[ThreadNum] + 1 : 0), NRows = 0, NCols = 0;
    Potential *Pot = (ThreadNum >= 0 ? Pots[ThreadNum] : OPot);
    QString **PotData = Pot->getData(NRows, NCols), Res[NRows + 2];
    for (n=0; n < NRows; n++) Res[n] = PotData[n][1];
    Res[n] = QString::number(SFQS, 'g', 10);
    Res[n+1] = QString::number(FQS, 'g', 12);
    ResultTab->setRowData(r, Res);
    Destroy(PotData, NRows);
    if (ThreadNum >= 0) ThreadControl::CalcFinished(ThreadNum);
    else 
    {
        Pot->getSFQSU(SFQSU, mJ, mv);
        NumItE->setText(QString::number(PotList.count()));
        NItChanged();
    }
}

void SFQSCalcControl::StartCalc(int p, int n)
{
    if (Pots[p] == 0)
    {
        Pots[p] = MW->CreatePotential(p);
        Pots[p]->show();
        connect(Pots[p], SIGNAL(FitTerminated(int)), this, SLOT(CalcTerminated(int)));
        connect(Pots[p], SIGNAL(SFQScalculated(int,double,double)), this, SLOT(CalcFinished(int,double,double)));
    }
    St->removePotential(Pots[p]);
    Pots[p]->readData(PDir + PotList[n]);
    if (!FDir.isEmpty())
    {
        MFD->readData(FDir + PotList[n].left(PotList[n].length() - 3) + "fdat");
        if (St->getFitData() != MFD) St->setMainFitData(MFD);
    }
    St->addPotential(Pots[p]);
    Pots[p]->calcFQS_SFQS(SFQSU);
}
