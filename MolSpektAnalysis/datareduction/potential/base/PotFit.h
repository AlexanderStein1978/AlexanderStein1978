//
// C++ Interface: PotFit
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POTFIT_H
#define POTFIT_H


#include <QThread>


class PotFit : public QThread
{
	Q_OBJECT
	
public:
	enum Task {PotentialFit, MonteCarloSimulation, calcFQS};
	
    PotFit();
	void FitAnaPot(bool robustWeighting = false, bool UseSVD = true, double threshhold = 1e-10,
				   bool useLevenbergMarquardt = true, int MaxIt = -1, bool adjustTEWeightFact = false, double hCi_cP = 0.0);
	void setWorker(PotWorker *Worker);
	
	inline void FitSplinePot(int nMaxIt = 1, double threshFakt = 0.9999)
	{
		ToDo = PotentialFit;
		FitPot = SplinePotential;
		MaxIt = nMaxIt;
		threshhold = threshFakt;
		start();
	}
	
	inline void CalcFQS(double nSFQSRad = 0.0)
	{
		ToDo = calcFQS;
		SFQSRad = nSFQSRad;
		start();
	}
	
	inline void restartFit()
	{
		restart = true;
		start();
	}
	
	inline void RunMonteCarloIteration(double *MCSTab, double UncFact = 1.0)
	{
		ToDo = MonteCarloSimulation;
		MCST = MCSTab;
		MCSUncFact = UncFact;
		start();
	}
	
	inline PotFit::Task getTask()
	{
		return ToDo;
	}
	
protected:
	void run();
	
private:
	
	Task ToDo;
	bool robustWeighting, UseSVD, adjustTEWeightFact, useLevenbergMarquardt, restart;
	int MaxIt;
	double threshhold, hCi_cP, *MCST, SFQSRad, MCSUncFact;
	PotWorker *Worker;
	PotentialType FitPot;
};

#endif
