//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FITDATA_BASEDATA_H
#define FITDATA_BASEDATA_H


#include "basedata.h"


struct FitDataBaseData : BaseData
{
	QString vs, source, secondState;
	double devR = 0.0f;
};

#endif
