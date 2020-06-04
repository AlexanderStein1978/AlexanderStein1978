//
// C++ Interface: Point
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POINT_H
#define POINT_H


struct Point
{
    double x, y, sig, unc, obs;
	int SplineNum, MarkerNum;
};

#endif
