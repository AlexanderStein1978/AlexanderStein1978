//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
