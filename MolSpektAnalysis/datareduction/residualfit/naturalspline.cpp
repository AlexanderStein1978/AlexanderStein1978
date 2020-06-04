//
// C++ Implementation: NaturalSpline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include <QList>
#include <QTextStream>
#include <QStringList>
#include <QMessageBox>

#include "naturalspline.h"
#include "molecule.h"
#include "elstate.h"
#include "utils.h"
#include "fit.h"


typedef QList<SplinePoint*> SplinePointArrayList;


NaturalSpline::NaturalSpline()
{
	L = 0;
	points = 0;
	numPoints = 0;
}

NaturalSpline::NaturalSpline(SplinePoint* Npoints, int N)
{
	L = 0;
	points = Npoints;
	numPoints = N;
	CalcLMatrix();
}

NaturalSpline::~NaturalSpline()
{
	if (points != 0) delete[] points;
	if (L != 0) Destroy(L, numPoints);
}

void NaturalSpline::CalcLMatrix()
{
	int N = numPoints, n, i, st = numPoints -1;
	if (N == 0) return;
	if (L==0) L = Create(N, N);
	for (n=0; n<N; n++) for (i=0; i<N; i++) L[n][i] = 0.0;
	points[0].yss = points[st].yss = 0.0;
	if (N<3) return;
	double B[st], F, dx[st];
	for (i=0; i < st; i++) dx[i] = points[i+1].x - points[i].x;
	B[1] = 1.0 / (points[2].x - points[0].x);
	for (n=1; n < st; n++) 
	{
		L[n][n-1] = (n>1 ? L[n-1][n] : 3.0 / dx[0]);
		L[n][n+1] = 3.0 / dx[n];
		L[n][n] = - L[n][n-1] - L[n][n+1];
	}
	for (n=2; n < st; n++)
	{
		F = 0.5 * dx[n-1] * B[n-1];
		B[n] = 1.0 / (points[n+1].x - points[n-1].x - F * 0.5 * dx[n-1]);
		for (i=0; i<=n; i++) L[n][i] -= F * L[n-1][i];
	}
	for (i=0; i<N; i++) L[st-1][i] *= B[st-1];
		//		(L[st][i] - 0.5 * dx[st + 1] * points[st+2].yss) * B[st - 1];
	for (n = st - 2; n>0; n--) for (i=0; i<N; i++) 
			L[n][i] = (L[n][i] - 0.5 * dx[n] * L[n+1][i]) * B[n];
	for (n=1; n < st; n++) for (i=0, points[n].yss = 0.0; i<N; i++) points[n].yss += L[n][i] * points[i].y;
}

double* NaturalSpline::getCurve(double xmin, double xmax, int NPoints)
{
	int n, i, x, N = numPoints, st = numPoints - 1;
	double *y = new double[NPoints];
	if (N<2)
	{
		if (N==0) for (n=0; n < NPoints; n++) y[n] = 0.0;
		else for (n=0; n < NPoints; n++) y[n] = points[0].y;
		return y;
	}
	double dx[st], A, B, h = (xmax - xmin) / double(NPoints - 1), es = 1.0 / 6.0;
	for (n=0; n < st; n++) dx[n] = points[n+1].x - points[n].x;
	for (x = xmin, n=0; x < points[0].x; x+=h, n++) y[n] = 0.0;
	for (i=0; n < NPoints && x <= points[st].x; n++, x+=h)
	{
		while(points[i].x < x) i++;
		if (i>0)
		{
			B = 1.0 - (A = (points[i].x - x) / dx[i-1]);
			y[n] = A * points[i-1].y + B * points[i].y + es * (dx[i-1] * dx[i-1]) * ((A*A*A-A) * points[i-1].yss 
				 + (B*B*B-B) * points[i].yss);
		}
		else y[n] = points[0].y;
	}
	while (n < NPoints) y[n++] = 0.0;
	if (xmax == points[st].x) y[NPoints - 1] = points[st].y; 
	return y;
}

double* NaturalSpline::getDerivation(double Rmin, double Rmax, int NPoints)
{
	if (points == 0) return 0;
	int n, m;
	double r, *R = new double[NPoints], h = (Rmax - Rmin) / (NPoints - 1), xd = 0.0, A, B, eds = 0.1666666666667, T = 0.0, F = 0.0;
	for (n=0, r = Rmin; n < NPoints && r < points[0].x; n++, r+=h) R[n] = 0.0;
	for (m=-1; n < NPoints && r < points[numPoints - 1].x; n++, r+=h)
	{
		if (r > points[m+1].x)
		{
			while (r > points[m+1].x) m++;
			xd = 1.0 / (points[m+1].x - points[m].x);
			T = (points[m+1].y - points[m].y) * xd;
			F = eds * (points[m+1].x - points[m].x);
		}
		B = (r - points[m].x) * xd;
		A = 1.0 - B;
		R[n] = T + F * ((1.0 - 3.0 * A * A) * points[m].yss + (3.0 * B * B - 1.0) * points[m+1].yss);
	}
	while (n < NPoints) R[n] = 0.0;
	return R;
}

double NaturalSpline::getDerivative(double R)
{
	if (points == 0) return 0.0;
	if (R < points[0].x || R > points[numPoints - 1].x) return 0.0;
	int m;
	double xd, A, B, eds = 0.1666666666667, T, F;
	for (m=0; R > points[m+1].x; m++) ;
	xd = 1.0 / (points[m+1].x - points[m].x);
	T = (points[m+1].y - points[m].y) * xd;
	F = eds * (points[m+1].x - points[m].x);
	B = (R - points[m].x) * xd;
	A = 1.0 - B;
	return T + F * ((1.0 - 3.0 * A * A) * points[m].yss + (3.0 * B * B - 1.0) * points[m+1].yss);
}

int NaturalSpline::getNumPoints()
{
	return numPoints;
}

double NaturalSpline::getPoint(double x)
{
	if (numPoints > 0 ? x < points[0].x || x > points[numPoints - 1].x : true) return 0.0;
	int i;
	for (i=0; i < numPoints; i++) if (x == points[i].x) return points[i].y;
	for (i=0; points[i].x < x; i++) ;
	double dx = points[i].x - points[i-1].x;
	double A = (points[i].x - x) / dx;
	double B = 1.0 - A;
	return A * points[i-1].y + B * points[i].y 
	         + dx * dx / 6.0 * ((A*A*A-A) * points[i-1].yss + (B*B*B-B) * points[i].yss); 
}

double** NaturalSpline::getS(double xmin, double xmax, int NPoints)
{
	int n, s, i, x, N = numPoints, st = numPoints - 1;
	if (N==0) return 0;
	double **S = Create(N, NPoints);
	if (N<2)
	{
		for (n=0; n < NPoints; n++) S[0][n] = 0.0;
		return S;
	}
	double dx[st], dxq, A, B, C, D, h = (xmax - xmin) / double(NPoints - 1), es = 1.0 / 6.0;
	for (n=0; n < st; n++) dx[n] = points[n+1].x - points[n].x;
	for (x = xmin, n=0; x < points[0].x; x+=h, n++) for (s=0; s<N; s++) S[s][n] = 0.0;
	for (i=0; n < NPoints && x <= points[st].x; n++, x+=h)
	{
		while(points[i].x < x) i++;
		if (i>0)
		{
			B = 1.0 - (A = (points[i].x - x) / dx[i-1]);
			C = (dxq =  es * (dx[i-1] * dx[i-1])) * (A*A*A-A);
			D = dxq * (B*B*B-B); 
			for (s=0; s<N; s++) S[s][n] = C * L[i-1][s] + D * L[i][s]; 
			S[i-1][n] += A;
			S[i][n] += B;
		}
		else for (s=1, S[0][n] = 1.0; s<N; s++) S[s][n] = 0.0;
	}
	for ( ; n < NPoints; n++) for (s=0; s<N; s++) S[s][n] = 0.0;
	if (xmax == points[st].x) for (s=0, S[st][NPoints - 1] = 1.0; s < st; s++) S[s][NPoints - 1] = 0.0; 
	return S;
}

void NaturalSpline::setPoints(SplinePoint* Npoints, int N)
{
	if (points != 0) delete[] points;
	if (L!=0) 
	{
		Destroy(L, numPoints);
		L=0;
	}
	points = Npoints;
	numPoints = N;
	CalcLMatrix();
}

void NaturalSpline::setyValues(double* yi)
{
	int n, i;
	for (n=0; n < numPoints; n++) points[n].y = yi[n];
	for (n=1; n < numPoints - 1; n++) for (i=0, points[n].yss = 0.0; i < numPoints - 1; i++) 
		points[n].yss += L[n][i] * points[i].y;
}
