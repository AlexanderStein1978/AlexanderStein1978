//
// C++ Implementation: ImprovePotSeriesControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "ImprovePotSeriesControl.h"


ImprovePotSeriesControl::ImprovePotSeriesControl(MainWindow* MW, ElState* St, Potential* OPot, FitData* FDat, QString potDir, 
                                                 double Threshold, int NumParFits) : ThreadControl(MW, St, 0, NumParFits)
{
    QDir Dir(PotDir = potDir);
    int NRows = 0, NCols = 0, c;
    threshold = Threshold;
    if (PotDir.right(1) != "\\" && PotDir.right(1) != "/") PotDir += DIRSEP;
    PotList = Dir.entryList(QStringList("*.pot"), QDir::Files);
    if (PotList.count() == 0)
    {
        CalcFinished(0);
        return;
    }
    if (OPot != 0) OPot->getFixedCoefficients(NFC, FCRows);
    else
    {
        OPot = Pots[0] = MW->CreatePotential(0);
        OPot->readData(PotDir + PotList[0]);
        Pots[0]->show();
        connect(Pots[0], SIGNAL(FitTerminated(int)), this, SLOT(CalcTerminated(int)));
        connect(Pots[0], SIGNAL(FitFinished(int, double)), this, SLOT(FitFinished(int,double)));
    }
    UseLeveMarq = OPot->getPotType() != SplinePotential;
    QString **PotData = OPot->getData(NRows, NCols);
    QStringList HeaderLabels;
    setWindowTitle("Improve Potential Series Control");
    AdTabRows = 0;
    ResultTab->setWindowTitle("Results of Potential Series Improvement");
    for (c=0; c < NRows; c++) HeaderLabels << PotData[c][0];
    HeaderLabels << "FQS" << "FQS_Bad" << "NBad" << "NBadPAL" << "Sigma";
    ResultTab->setTabDimensions(PotList.count() + 1, NumResultTabColumns = NRows + 5);
    ResultTab->setHorizontalHeader(HeaderLabels);
    Destroy(PotData, NRows);
    if (St->getFitData() != FDat) St->setMainFitData(FDat);
    NumItE->setText(QString::number(PotList.count()));
    NItChanged();
}

void ImprovePotSeriesControl::FitFinished(int ThreadNum, double FQS)
{
    int n, r = PotI[ThreadNum], NRows = 0, NCols = 0;
    Potential *Pot = Pots[ThreadNum];
    QString **PotData = Pot->getData(NRows, NCols), Res[NRows + 5], FileName = Pot->getFileName();
    FitResult FRes = Pot->getFitResult();
    for (n=0; n < NRows; n++) Res[n] = PotData[n][1];
    Res[n] = QString::number(FQS, 'g', 10);
    Res[n+1] = QString::number(FRes.FQS_Bad, 'g', 10);
    Res[n+2] = QString::number(FRes.NBad);
    Res[n+3] = QString::number(FRes.NBadPAL);
    Res[n+4] = QString::number(FRes.Sigma, 'g', 10);
    ResultTab->setRowData(r, Res);
    Destroy(PotData, NRows);
    Pot->writeData(getWriteFileName(FileName));
    CalcFinished(ThreadNum);
}

QString ImprovePotSeriesControl::getWriteFileName(QString CurrentFileName)
{
    return CurrentFileName.left(CurrentFileName.length() - 4) + "_i.pot";
}

void ImprovePotSeriesControl::StartCalc(int p, int n)
{
    PreparePot(p, n);
    Pots[p]->FitAnaPot(false, false, false, threshold, false, false, true, -1, false, 0.0, true);
}

void ImprovePotSeriesControl::PreparePot(int p, int n)
{
    if (Pots[p] == 0)
    {
        Pots[p] = MW->CreatePotential(p);
        Pots[p]->show();
        connect(Pots[p], SIGNAL(FitTerminated(int)), this, SLOT(CalcTerminated(int)));
        connect(Pots[p], SIGNAL(FitFinished(int,double)), this, SLOT(FitFinished(int,double)));
    }
    St->removePotential(Pots[p]);
    Pots[p]->readData(PotDir + PotList[n]);
    if (NFC > 0) Pots[p]->VaryCoefficients(false, NFC, FCRows);
    St->addPotential(Pots[p]);
}
