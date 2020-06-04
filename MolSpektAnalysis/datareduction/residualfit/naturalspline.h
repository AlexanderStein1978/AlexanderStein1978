//
// C++ Interface: NaturalSpline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef NATURALSPLINE_H
#define NATURALSPLINE_H


#include <vector>
#include <QObject>

#include "fitobject.h"

template<class T> class QList;
class QTextStream;
class QString;

class ElState;


class NaturalSpline
{
	public:
		NaturalSpline();
		NaturalSpline(SplinePoint *points, int N);
		~NaturalSpline();
		void setPoints(SplinePoint *points, int N);
		void setyValues(double *yi);
		int getNumPoints();
		double getPoint(double R);
		double getDerivative(double R);
		double *getDerivation(double Rmin, double Rmax, int NPoints);
		double *getCurve(double Rmin, double Rmax, int NPoints);
		double **getS(double Rmin, double Rmax, int NPoints);
		
	private:
		void CalcLMatrix();
		
		int numPoints;
		double **L;
		SplinePoint *points;
};

#endif // NATURALSPLINE_H
