//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
