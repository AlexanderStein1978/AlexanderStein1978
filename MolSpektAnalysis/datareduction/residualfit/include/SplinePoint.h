//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SPLINEPOINT_H
#define SPLINEPOINT_H


struct SplinePoint
{
	double x, y, yss;
	bool variable;
	
	inline SplinePoint& operator+=(SplinePoint &other)
	{
		x += other.x;
		y += other.y;
		yss += other.yss;
		return *this;
	}
	
	inline SplinePoint& operator/=(double div)
	{
		x /= div;
		y /= div;
		yss /= div;
		return *this;
	}
	
	inline bool operator!=(SplinePoint &other)
	{
		return x != other.x || y != other.y || yss != other.yss;
	}
};

#endif
