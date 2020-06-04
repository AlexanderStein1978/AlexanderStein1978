//
// C++ Implementation: SplinePot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "splinepot.h"


SplinePot::SplinePot(PotFit* Fit): PotWorker(Fit, SplinePotential)
{

}

SplinePot::SplinePot(const SplinePot& C): PotWorker(C)
{

}
