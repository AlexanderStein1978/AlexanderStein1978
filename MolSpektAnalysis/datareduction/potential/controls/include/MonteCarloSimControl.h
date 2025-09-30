//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MONTECARLOSIMCONTROL_H
#define MONTECARLOSIMCONTROL_H


#include "ThreadControl.h"


class MonteCarloSimControl : public ThreadControl
{
	Q_OBJECT
	
public:
    MonteCarloSimControl(MainWindow* MW, Potential *Pot, QString PotDir, int NumIterations, int NumParFits, double UncFact, int NumWFPoints);
    ~MonteCarloSimControl();
	
public slots:
	void FitFinished(int ThreadNum, double FQS);
	
protected:
	void StartCalc(int p, int n);
	
private:
	int NumPDRows, NumPDColumns;
	double *MCSTab, UncF;
	FitData *MCSFD;
};

#endif
