//
// C++ Implementation: ConnectSpline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "ConnectSpline.h"


ConnectSpline::ConnectSpline(const Spline * const i_leftSpline, const Spline * const i_rightSpline)
    : Spline(2, new SplinePoint[2], i_leftSpline == 0 || i_rightSpline == 0), m_leftSpline(i_leftSpline), m_rightSpline(i_rightSpline)
{
    Interpolate();
}

void ConnectSpline::Interpolate()
{
    if (m_leftSpline != 0)
    {
        points[0].x = m_leftSpline->getxN();
        points[0].y = m_leftSpline->gety(points[0].x);
    }
    if (m_rightSpline != 0)
    {
        points[1].x = m_rightSpline->getx0();
        points[1].y = m_rightSpline->gety(points[1].x);
    }
    if (m_leftSpline != 0 && m_rightSpline != 0)
    {
        const double ddx = 1.0 / (points[1].x - points[0].x), dyBlock = 6.0 * (points[1].y - points[0].y) * ddx * ddx;
        const double ysl = m_leftSpline->getys(points[0].x), ysr = m_rightSpline->getys(points[1].x);
        points[1].yss = 2.0 * (ysl + 2.0 * ysr) * ddx - dyBlock;
        points[0].yss = 6.0 * ysr * ddx - dyBlock - 2.0 * points[1].yss;
    }
    else
    {
        points[0].yss = points[1].yss = 0.0;
        if (m_rightSpline != 0)
        {
            const double ysr = m_rightSpline->getys(points[1].x);
            points[0].x = 0.0;
            points[0].y = points[1].y - ysr * points[1].x;
        }
        else if (m_leftSpline != 0)
        {
            const double ysl = m_leftSpline->getys(points[0].x);
            points[1].x = 1e3;
            points[1].y = points[0].y + ysl * (points[1].x - points[0].x);
        }
        else points[0].x = points[1].x = points[0].y = points[1].y = 0.0;
    }
}
