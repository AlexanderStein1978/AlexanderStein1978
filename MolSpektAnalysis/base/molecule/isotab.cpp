//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "isotab.h"

#include <QPixmap>

#include <cmath>


IsoTab::IsoTab()
{
}

IsoTab::~IsoTab()
{
	//printf("IsoTab::~IsoTab()\n");
	if (mNumIso1 != nullptr) delete[] mNumIso1;
	if (mNumIso2 != nullptr) delete[] mNumIso2;
	if (relNA != nullptr) delete[] relNA;
	if (redMass != nullptr) delete[] redMass;
	if (relRedMass != nullptr) delete[] relRedMass;
	if (rootRRM != nullptr) delete[] rootRRM;
	if (texName != nullptr) delete[] texName;
	if (JStep != nullptr) delete[] JStep;
	if (mIso1 != nullptr) delete[] mIso1;
	if (mIso2 != nullptr) delete[] mIso2;
	if (chSymb1 != nullptr) delete chSymb1;
	if (chSymb2 != nullptr) delete chSymb2;
	if (IsoImage != nullptr) delete[] IsoImage;
	if (Isoname != nullptr) delete[] Isoname;
	//printf("Ende ~IsoTab\n");
}

void IsoTab::createArrays()
{
	IsoImage = new QPixmap[numIso];
	Isoname = new QString[numIso];
	JStep = new int[numIso];
	mIso1 = new double[numIso];
	mIso2 = new double[numIso];
	mNumIso1 = new int[numIso];
	mNumIso2 = new int[numIso];
	redMass = new double[numIso];
	relNA = new double[numIso];
	relRedMass = new double[numIso];
	rootRRM = new double[numIso];
	texName = new QString[numIso];
	chSymb1 = new QString;
	chSymb2 = new QString;
}

int IsoTab::getIsoIndex(int mI1, int mI2)
{
    for (int n=0; n < numIso; ++n) if ((mNumIso1[n] == mI1 && mNumIso2[n] == mI2) || (mNumIso1[n] == mI2 && mNumIso2[n] == mI1)) return n;
    return -1;
}
