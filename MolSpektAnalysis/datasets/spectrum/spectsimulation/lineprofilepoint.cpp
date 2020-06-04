//
// C++ Implementation: LineProfilePoint
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "lineprofilepoint.h"


LineProfilePoint::LineProfilePoint(const double i_Energy, const double i_Intensity)
    : Energy(i_Energy), Intensity(i_Intensity)
{
}

LineProfilePoint::LineProfilePoint()
    : Energy(0.0), Intensity(0.0)
{
}
