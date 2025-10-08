//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "fitresult.h"

FitResult::FitResult()
{
	initialize();
}

FitResult::FitResult(const FitResult& toCopy)
 : NBad(toCopy.NBad), NBadPAL(toCopy.NBadPAL), Index(toCopy.Index), ProcessNum(toCopy.ProcessNum)
 , Sigma(toCopy.Sigma), FQS_Bad(toCopy.FQS_Bad), ParValue(toCopy.ParValue), LRelParDev(toCopy.LRelParDev)
 , isStartPot(toCopy.isStartPot)
 {
 }
 
 FitResult::~FitResult()
 {
 }

void FitResult::initialize()
{
	FQS_Bad = Sigma = -1.0;
	Index = NBad = NBadPAL = ProcessNum = -1;
	ParValue = LRelParDev = 0.0;
	isStartPot = false;
}

bool FitResult::isValid()
{
	if (FQS_Bad < 0.0 || Sigma <= 0.0 || NBad < 0 || NBadPAL < 0) return false;
	return true;
}

bool FitResult::better(FitResult *Result2)
{
	if (LRelParDev < 10.0 && Result2->LRelParDev >= 10.0) return true;
	if (LRelParDev >= 10.0 && Result2->LRelParDev < 10.0) return false;
	if (NBadPAL < Result2->NBadPAL) return true;
	if (NBadPAL > Result2->NBadPAL) return false;
	if (NBad < Result2->NBad) return true;
	if (NBad > Result2->NBad) return false;
	if (FQS_Bad < Result2->FQS_Bad) return true;
	if (FQS_Bad > Result2->FQS_Bad) return false;
	if (Sigma < Result2->Sigma) return true;
	return false;
}

bool FitResult::betterEqual(FitResult *Result2)
{
	if (NBadPAL < Result2->NBadPAL) return true;
	if (NBadPAL > Result2->NBadPAL) return false;
	if (NBad < Result2->NBad) return true;
	if (NBad > Result2->NBad) return false;
	if (FQS_Bad < Result2->FQS_Bad) return true;
	if (FQS_Bad > Result2->FQS_Bad) return false;
	if (Sigma <= Result2->Sigma) return true;
	return false;
}

bool FitResult::equal(FitResult* Result2)
{
	if (NBadPAL == Result2->NBadPAL && NBad == Result2->NBad && FQS_Bad == Result2->FQS_Bad && Sigma == Result2->Sigma) 
		return true;
	return false;
}
