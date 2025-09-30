//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "CreateAnaPotSeriesControl.h"
#include "tablewindow.h"
#include "elstate.h"
#include "molecule.h"
#include "fitdata.h"
#include "potential.h"

#include <QTextStream>
#include <QFile>
#include <QDir>



CreateAnaPotSeriesControl::CreateAnaPotSeriesControl(MainWindow* MW, ElState* St, Potential* OPot, FitData* fDat, QString SPDir, QString FDDir, QString APDir, 
                                                     double Threshold, int NumParFits, int maxFitIterations, int nStart, bool ImproveAnaPots, bool useSvd, bool useLeveMarq) 
                                                            : ImprovePotSeriesControl(MW, St, OPot, fDat, SPDir, Threshold, 0)
{
    setWindowTitle("Create AnaPot Series Control");
    ResultTab->setWindowTitle("Results of creating of AnaPot series");
    FitDataDir = FDDir;
    ResultPotDir = APDir;
    MaxFitIterations = maxFitIterations;
    NStart = nStart;
    improveAnaPots = ImproveAnaPots;
    UseSvd = useSvd;
    UseLeveMarq = useLeveMarq;
    FDat = fDat;
    PotN = St->getNumPotentials();
    FDN = St->getNumFitDataSets();
    MFD = St->getFitDataFileName();
    if (improveAnaPots) SPDir = APDir;
    FDat = (FDDir.isEmpty() ? 0 : MW->CreateFitData());
    if (APDir.right(1) != "/" && APDir.right(1) != "\\") APDir += '/';
    if (FDat != 0 && FDDir.right(1) != "/" && FDDir.right(1) != "\\") FDDir += '/';
    LogFile = new QFile(APDir + "CreateAnaPots.log");
    LogFile->open(QIODevice::WriteOnly);
    LogStream = new QTextStream(LogFile);
    *LogStream << "CreateAnaPotSeriesFromMCSSplinePotSeries(\"" << St->getMolecule()->getFileName() << "\", \"" << St->getName() << "\", \"" 
               << FDDir << "\", \"" << SPDir << "\", \"" << APDir << (improveAnaPots ? "\", true, " : "\", false, ")
               << (UseSvd ? "true, " : "false, ") << (UseLeveMarq ? "true, " : "false , ") << QString::number(MaxFitIterations) 
               << ", " << QString::number(Threshold, 'g', 5) << ")\n";
    LogStream->flush();
    NumParFE->setText(QString::number(NumParFits));
    NParFChanged();
}

CreateAnaPotSeriesControl::~CreateAnaPotSeriesControl()
{
    delete LogStream;
    delete LogFile;
}

void CreateAnaPotSeriesControl::AllCalculationsFinished()
{
    *LogStream << "Series finished!!!\n";
    St->removePotential(PotN);
    if (FDat != 0 && (UseLeveMarq || UseSvd))
    {
        St->setMainFitData(MFD);
        St->removeFitData(FDN);
    }
    ThreadControl::AllCalculationsFinished();
}

QString CreateAnaPotSeriesControl::getWriteFileName(QString CurrentFileName)
{
    return ImprovePotSeriesControl::getWriteFileName(CurrentFileName);
}

void CreateAnaPotSeriesControl::StartCalc(int p, int n)
{
    PreparePot(p, n);
    QDir SplinePotDir(PotDir);
    QString SFile, CFile = SplinePotDir.absoluteFilePath(PotList[n]);
    if (!improveAnaPots)
    {
        SFile = ResultPotDir + PotList[n].left(PotList[n].lastIndexOf('.')) + "AnaPot.pot";
        QFile File(SFile);
        if (File.exists()) 
        {
            CalcFinished(p);
            return;
        }
    }
    if (improveAnaPots) SFile = CFile;
    printf("Potential%d: %s\n", n, CFile.toLatin1().data());
    if (FDat != 0 && (UseLeveMarq || UseSvd))
    {
        if (improveAnaPots)
        {
            CFile = FitDataDir + PotList[n].left(PotList[n].indexOf("AnaPot.pot")) + ".fdat";
            QFile File(CFile);
            if (!File.exists()) CFile = FitDataDir + PotList[n].left(PotList[n].lastIndexOf('.') + 1) + "fdat";
        }
        else CFile = FitDataDir + PotList[n].left(PotList[n].lastIndexOf('.') + 1) + "fdat";
        if (!FDat->readData(CFile))
        {
            *LogStream << "Error reading file \"" + CFile + "\"!";
            CalcFinished(p);
            return;
        }
        if (St->getFitData() != FDat) St->setMainFitData(FDat);
    }
    if (!improveAnaPots)
    {
        Pots[p]->setFileName(SFile);
        Pots[p]->setName(Pots[p]->getName() + "AnaPot");
        Pots[p]->setSource("SplinePotential from " + Pots[p]->getSource());
        Pots[p]->FitAnaPot(false, false, UseSvd, threshold, true, false, UseLeveMarq);
    }
    else if (UseLeveMarq || UseSvd) Pots[p]->FitAnaPot(false, false, UseSvd, threshold, false, false, UseLeveMarq, MaxFitIterations);
}

void CreateAnaPotSeriesControl::CalcFinished(int ThreadNum)
{
    QString SFile = Pots[ThreadNum]->getFileName();
    int n = ItFinished;
    ThreadControl::CalcFinished(ThreadNum);
    *LogStream << "Run" << QString::number(n) << " with file \"" << SFile << "\" finished!\n";
    LogStream->flush();
}
