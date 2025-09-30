//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINEPROFILEPOINT_H
#define LINEPROFILEPOINT_H


struct LineProfilePoint
{
    LineProfilePoint();
    LineProfilePoint(const double i_Energy, const double i_Intensity);

    double Energy;
    double Intensity;
};

#endif
