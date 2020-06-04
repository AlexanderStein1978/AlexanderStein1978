//
// C++ Interface: SFQSCalcControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SFQSCALCCONTROL_H
#define SFQSCALCCONTROL_H


#include "ThreadControl.h"


class SFQSCalcControl : public ThreadControl
{
	Q_OBJECT
	
public:
	SFQSCalcControl(MainWindow *MW, ElState *St, Potential *OPot, QString PotDir, QString FDatDir, double SFQSRad, int NumParFits); 
    virtual ~SFQSCalcControl();
	
public slots:
	void CalcFinished(int ThreadNum, double SFQS, double FQS);
	
protected:
	void StartCalc(int p, int n);
	
private:
	QString PDir, FDir;
	QStringList PotList;
	double Rad;
	bool ****SFQSU;
	int *mJ, **mv;
	Potential *OPot;
};

#endif
