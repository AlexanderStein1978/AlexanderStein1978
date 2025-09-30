//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
