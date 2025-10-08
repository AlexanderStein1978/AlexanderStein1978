//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CONNECTSPLINE_H
#define CONNECTSPLINE_H


#include "Spline.h"


class ConnectSpline : private Spline
{
public:
    ConnectSpline(const Spline* const i_leftSpline, const Spline* const i_rightSpline);
    void Interpolate();

    inline double GetValue(const double JValue) const
    {
        return gety(JValue);
    }

    inline double GetJMin()
    {
        return getx0();
    }

    inline double GetJMax()
    {
        return getxN();
    }

private:
    const Spline *const m_leftSpline, *const m_rightSpline;
};

#endif
