//
// C++ Implementation: IsoTab
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "isotab.h"


IsoTab::IsoTab(int NumIso)
{
	numIso = NumIso;
	mNumIso1 = new int[numIso];
	mNumIso2 = new int[numIso];
	relNA = new double[numIso];
	redMass = new double[numIso];
	relRedMass = new double[numIso];
	rootRRM = new double[numIso];
	texName = new QString[numIso];
	JStep = new int[numIso];
	mIso1 = new double[numIso];
	mIso2 = new double[numIso];
	chSymb1 = new QString;
	chSymb2 = new QString;
}

IsoTab::~IsoTab()
{
	//printf("IsoTab::~IsoTab()\n");
	delete[] mNumIso1;
	delete[] mNumIso2;
	delete[] relNA;
	delete[] redMass;
	delete[] relRedMass;
	delete[] rootRRM;
	delete[] texName;
	delete[] JStep;
	delete[] mIso1;
	delete[] mIso2;
	delete chSymb1;
	delete chSymb2;
	//printf("Ende ~IsoTab\n");
}

QString IsoTab::getIsoName(int N)
{
	if (N < 0 || N >= numIso) return "";
	if (mNumIso1[N] == mNumIso2[N] && *chSymb1 == *chSymb2) return QString::number(mNumIso1[N]) + *chSymb1 + "2";
	return QString::number(mNumIso1[N]) + *chSymb1 + QString::number(mNumIso2[N]) + *chSymb2;
}
