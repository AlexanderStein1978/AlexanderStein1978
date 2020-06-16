//
// C++ Interface: SplinePot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SPLINEPOT_H
#define SPLINEPOT_H


#include "potworker.h"


class SplinePot : public PotWorker
{
public:
	SplinePot(PotFit *Fit);
	SplinePot(const SplinePot &C);
};

#endif
