//
// C++ Interface: ImprovePotSeriesControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef IMPROVEPOTSERIESCONTROL_H
#define IMPROVEPOTSERIESCONTROL_H


#include "ThreadControl.h"


class ImprovePotSeriesControl : public ThreadControl
{
	Q_OBJECT
	
public:
	ImprovePotSeriesControl(MainWindow *MW, ElState* St, Potential *OPot, FitData *FDat, QString PotDir, double Threshold, 
							int NumParFits);
	
public slots:
	void FitFinished(int ThreadNum, double FQS);
	
protected:
	void StartCalc(int p, int n);
	void PreparePot(int p, int n);
	virtual QString getWriteFileName(QString CurrentFileName);
	
	QString PotDir;
	QStringList PotList;
	double threshold;
	bool UseLeveMarq;
};

#endif
