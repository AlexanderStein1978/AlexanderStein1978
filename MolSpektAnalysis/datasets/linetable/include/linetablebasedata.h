//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LINETABLE_BASEDATA_H
#define LINETABLE_BASEDATA_H


#include "basedata.h"


struct LineTableBaseData : BaseData
{
    int vs = -1, vss = -1, Jss = -1, F = -1;
	double waveNumber = 0.0, FCF = 0.0, upperEnergy = 0.0, averageUpperEnergy = 0.0, diffToAverageEnergy = 0.0;
	double energyDiffToNextJ = 0.0, calculatedEnergy = 0.0, diffToCalculatedEnergy = 0.0;
	QString Comment;
	const QPixmap* IsoIcon = nullptr;
};

#endif
