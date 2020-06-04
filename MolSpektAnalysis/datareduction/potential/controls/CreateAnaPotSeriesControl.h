//
// C++ Interface: CreateAnaPotSeriesControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CREATEANAPOTSERIESCONTROL_H
#define CREATEANAPOTSERIESCONTROL_H


#include "ImprovePotSeriesControl.h"


class CreateAnaPotSeriesControl : public ImprovePotSeriesControl
{
public:
    CreateAnaPotSeriesControl(MainWindow* MW, ElState* St, Potential* OPot, FitData* FDat, QString SplinePotDir, QString FitDataDir, QString ResultPotDir, 
							  double Threshold, int NumParFits, int MaxFitIterations, int NStart, bool improveAnaPots, bool UseSvd, bool UseLeveMarq);
    virtual ~CreateAnaPotSeriesControl();
	
protected:
	void CalcFinished(int ThreadNum);
	void StartCalc(int p, int n);
	void AllCalculationsFinished();
	QString getWriteFileName(QString CurrentFileName);
	
private:
	QString FitDataDir, ResultPotDir, MFD;
	int MaxFitIterations, NStart, PotN, FDN;
	QTextStream *LogStream;
	QFile *LogFile;
	bool improveAnaPots, UseSvd, UseLeveMarq;
	FitData *FDat;
};

#endif
