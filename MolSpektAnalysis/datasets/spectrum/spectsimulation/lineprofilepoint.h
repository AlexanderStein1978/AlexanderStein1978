//
// C++ Interface: LineProfilePoint
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
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
