//
// C++ Interface: BandConst
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef BANDCONST_H
#define BANDCONST_H


struct BandConst
{
    BandConst() : value(0.0), stdDev(0.0)
    {
    }

    BandConst(const double i_value, const double i_stdDev) : value(i_value), stdDev(i_stdDev)
    {
    }

    double value, stdDev;
};

#endif
