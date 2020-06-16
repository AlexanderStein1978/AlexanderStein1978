//
// C++ Interface: FitResult
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FITRESULT_H
#define FITRESULT_H


#include <QMetaType>


struct FitResult
{
	FitResult();
	FitResult(const FitResult& toCopy);
	virtual ~FitResult();
	void initialize();
	bool better(FitResult *Result2);
	bool equal(FitResult *Result2);
	bool betterEqual(FitResult *Result2);
	bool isValid();
	
	int NBad, NBadPAL, Index, ProcessNum;
	double Sigma, FQS_Bad, ParValue, LRelParDev;
	bool isStartPot;
};

Q_DECLARE_METATYPE(FitResult)



#endif
