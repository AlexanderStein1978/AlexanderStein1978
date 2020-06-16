//
// C++ Implementation: Spline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "Spline.h"
#include "point.h"
#include "utils.h"
#include "LocalPerturbation.h"

#include <QTextStream>

#include "fit.h"


Spline::Spline() : L(0), S(0), S0(0.0), SN(0.0), NDataPoints(0), NSplinePoints(0), m_error(SENoError), DataPoints(0), points(0), Natural(false)
{
}

Spline::Spline(int nSplinePoints, SplinePoint* SPoints, bool natural)
    : L(0), S(0), S0(0.0), SN(0.0), NDataPoints(0), NSplinePoints(nSplinePoints), m_error(SENoError),
      DataPoints(0), points(SPoints), Natural(natural)
{
}

Spline::Spline(const Spline &i_other) : L(0), S(0), S0(i_other.S0), SN(i_other.SN),
    NDataPoints(i_other.NDataPoints), NSplinePoints(i_other.NSplinePoints), m_error(i_other.m_error), DataPoints(0), points(0),
    Natural(i_other.Natural)
{
    copyFields(i_other);
}

Spline::~Spline()
{
    if (points != 0) delete[] points;
    if (NDataPoints > 0) delete[] DataPoints;
	deleteSL();
}

void Spline::calcS()
{
    if (NSplinePoints < 2) return;
    int n, m, s = (Natural ? 0 : 1), p = 0;
	double A, B, C, D, DR, d6 = 1.0 / 6, DR2;
	if (L == 0)
	{
		L = Create(NSplinePoints, NSplinePoints + 2*s);
		CalcLMatrix(L, points, 0.0, NSplinePoints, Natural);
	}
	S = Create(NSplinePoints + 2*s, NDataPoints);
	for (n = 0; n < NDataPoints; n++)
	{
		while (DataPoints[n].x > points[p+1].x && p < NSplinePoints - 2) p++;
		DR = points[p+1].x - points[p].x;
		DR2 = DR * DR * d6;
		A = (points[p+1].x - DataPoints[n].x) / DR;
		B = 1.0 - A;
		C = A * (A*A - 1.0) * DR2;
		D = B * (B*B - 1.0) * DR2;
		for (m=0; m < NSplinePoints + 2 * s; m++) S[m][n] = L[p][m] * C + L[p+1][m] * D;
		S[p+s][n] += A;
		S[p+s+1][n] += B;
	}
}

void Spline::calcFQSRes(double* Res, double& FQS)
{
	int n;
	for (n=0, FQS = 0.0; n < NDataPoints; n++)
	{
		Res[n] = DataPoints[n].y - gety(DataPoints[n].x);
		FQS += Res[n] * Res[n] / (DataPoints[n].sig * DataPoints[n].sig);
	}
	printf("FQS=%f\n", FQS);
}

void Spline::calcYss()
{
	int n, m, s = (Natural ? 0 : 1);
	if (Natural && NSplinePoints <= 2) for (n=0; n < NSplinePoints; n++) points[n].yss = 0.0;
	else for (n=0; n < NSplinePoints; n++)
	{
		if (!Natural) points[n].yss = S0 * L[n][0];
		else points[n].yss = 0.0;
		for (m=0; m < NSplinePoints; m++) points[n].yss += L[n][m+s] * points[m].y;
		if (!Natural) points[n].yss += SN * L[n][m+s];
	}
}

void Spline::copyData(const int N, const Point *const dpoints)
{
	int n;
	if (DataPoints != 0) delete[] DataPoints;
	DataPoints = new Point[N];
	for (n=0; n<N; n++) DataPoints[n] = dpoints[n];
	NDataPoints = N;
}

void Spline::copyFields(const Spline &i_other)
{
    int n, m, s = (i_other.Natural ? 0 : 1);
    if (i_other.L != 0)
    {
        L = Create(i_other.NSplinePoints, i_other.NSplinePoints + 2*s);
        for (n=0; n < i_other.NSplinePoints; ++n) for (m=0; m < i_other.NSplinePoints + 2*s; ++m)
            L[n][m] = i_other.L[n][m];
    }
    if (i_other.S != 0)
    {
        S = Create(i_other.NSplinePoints + 2*s, i_other.NDataPoints);
        for (n=0; n < i_other.NSplinePoints + 2*s; ++n) for (m=0; m < i_other.NDataPoints; ++m)
            S[n][m] = i_other.S[n][m];
    }
    if (i_other.DataPoints != 0)
    {
        DataPoints = new Point[NDataPoints];
        for (n=0; n < NDataPoints; ++n) DataPoints[n] = i_other.DataPoints[n];
    }
    if (i_other.points != 0)
    {
        points = new SplinePoint[NSplinePoints];
        for (n=0; n < NSplinePoints; ++n) points[n] = i_other.points[n];
    }
}

void Spline::deleteSL()
{
	if (L != 0) Destroy(L, NSplinePoints);
	if (S != 0) Destroy(S, NSplinePoints + (Natural ? 0 : 2));
	S = L = 0;
}

void Spline::FitTo(int N, Point* dpoints)
{
	int n = 0;
	if (N == 0)
	{
		for (n=0; n < NSplinePoints; n++) 
			points[n].x = points[n].y = points[n].yss = 0.0;
		return;
	}
	if (points[0].x == points[NSplinePoints - 1].x)
	{
		double Y = 0.0, WS = 0.0, W;
		for (n=0; n<N; n++) 
		{
			W = dpoints[n].sig * dpoints[n].sig;
			Y += W * dpoints[n].y;
			WS += W;
		}
		Y /= WS;
		for (n=0; n < NSplinePoints; n++)
		{
			points[n].y = Y;
			points[n].yss = 0.0;
		}
		copyData(N, dpoints);
		return;
	}
	if (N < NSplinePoints)
	{
		for (n=0; n < NSplinePoints; n++)
		{
			points[n].x = dpoints[0].x;
			points[n].y = dpoints[0].y;
			points[n].yss = 0.0;
		}
		copyData(N, dpoints);
		return;
	}
	if (N == NSplinePoints)
	{
		for (n=0; n<N; n++) 
		{
			points[n].x = dpoints[n].x;
			points[n].y = dpoints[n].y;
			points[n].yss = 0.0;
		}
		copyData(N, dpoints);
		return;
	}
    if (S != 0)
    {
        if (N == NDataPoints) while (n < N && dpoints[n].x == DataPoints[n].x) n++;
        if (N != NDataPoints || n < N)
        {
            Destroy(S, Natural ? NSplinePoints : NSplinePoints + 2);
            S = 0;
        }
    }
    copyData(N, dpoints);
    Refit();
}

void Spline::Refit()
{
    if (NSplinePoints < 2) return;
    if (Natural && NSplinePoints == 2)
    {
        double *x = new double[NDataPoints], *y = new double[NDataPoints], *unc = new double[NDataPoints];
        for (int n=0; n < NDataPoints; n++)
        {
            x[n] = DataPoints[n].x;
            y[n] = DataPoints[n].y;
            unc[n] = DataPoints[n].unc;
        }
        FitToStraitLine(NDataPoints, x, y, unc, points[0].x, points[1].x, points[0].y, points[1].y);
        points[0].yss = points[1].yss = 0.0;
        delete[] x;
        delete[] y;
        delete[] unc;
        return;
    }
    int n, m, NC = (Natural ? NSplinePoints : NSplinePoints + 2), s = (Natural ? 0 : 1);
    double C[NC], *Y = new double[NDataPoints], *Sig = new double[NDataPoints], FQS, cFQS;
    if (S == 0) calcS();
    double **EQS = Create(NDataPoints,  NC);
    for (n=0; n < NDataPoints; n++)
    {
        Y[n] = DataPoints[n].y;
        Sig[n] = DataPoints[n].sig;
        for (m=0; m < NC; m++) EQS[n][m] = S[m][n];
    }
    SvdFit(EQS, C, Y, Sig, NC, NDataPoints, 1e-12);
    if (!Natural) S0 = C[0];
    for (n=0; n < NSplinePoints; n++) points[n].y = C[n+s];
    if (!Natural) SN = C[NC - 1];
    calcYss();
    calcFQSRes(Y, FQS);
    cFQS = 0.9 * FQS;
    while (FQS - cFQS > 1e-10)
    {
        FQS = cFQS;
        SvdFit(EQS, C, Y, Sig, NC, NDataPoints, 1e-12);
        UpdateY(C);
        calcYss();
        calcFQSRes(Y, cFQS);
    }
    delete[] Y;
    delete[] Sig;
}

void Spline::JoinForDeperturbation(Spline *i_other, int JminPert, int JmaxPert)
{
    int n, m, m1, m2, NewNumSplinePoints, NewNumDataPoints;
    for (NewNumSplinePoints = 0; NewNumSplinePoints < NSplinePoints && points[NewNumSplinePoints].x < JminPert; ++NewNumSplinePoints) ;
    for (m2=0; m2 < i_other->NSplinePoints && i_other->points[m2].x <= JmaxPert; ++m2) ;
    m1 = NewNumSplinePoints;
    NewNumSplinePoints += i_other->NSplinePoints - m2;
    SplinePoint* newPoints = new SplinePoint[NewNumSplinePoints];
    for (n=0; n < m1; ++n) newPoints[n] = points[n];
    for (m=m2; m < i_other->NSplinePoints; ++m, ++n) newPoints[n] = i_other->points[m];
    Natural = Natural && i_other->Natural;
    delete[] points;
    NSplinePoints = NewNumSplinePoints;
    points = newPoints;
    for (m1=0; m1 < NDataPoints && DataPoints[m1].x < JminPert; ++m1) ;
    for (m2=0; m2 < i_other->NDataPoints && i_other->DataPoints[m2].x <= JmaxPert; ++m2) ;
    NewNumDataPoints = m1 + i_other->NDataPoints - m2;
    Point* newDataPoints = new Point[NewNumDataPoints];
    for (n=0; n < m1; ++n) newDataPoints[n] = DataPoints[n];
    for (m=m2; m < i_other->NDataPoints; ++m, ++n) newDataPoints[n] = i_other->DataPoints[m];
    for (std::vector<LocalPerturbation*>::const_iterator it = i_other->m_localPerturbations.begin(); it != i_other->m_localPerturbations.end(); ++it) m_localPerturbations.push_back(*it);
    FitTo(NewNumDataPoints, newDataPoints);
}

Spline& Spline::operator =(const Spline& i_right)
{
    if (this != &i_right)
    {
        L = 0;
        S = 0;
        S0 = i_right.S0;
        SN = i_right.SN;
        NDataPoints = i_right.NDataPoints;
        NSplinePoints = i_right.NSplinePoints;
        DataPoints = 0;
        points = 0;
        Natural = i_right.Natural;
        copyFields(i_right);
    }
    return *this;
}

void Spline::addLocalPerturbation(LocalPerturbation* i_perturbation)
{
    m_localPerturbations.push_back(i_perturbation);
    for (std::vector<LocalPerturbation*>::reverse_iterator it = m_localPerturbations.rbegin(); it + 1 != m_localPerturbations.rend() && (*it)->GetCenter() < (*(it+1))->GetCenter(); ++it)
        std::swap(*it, *(it+1));
}

bool Spline::containsLocalPerturbation(const LocalPerturbation * const i_perturbation) const
{
    for (std::vector<LocalPerturbation*>::const_iterator it = m_localPerturbations.begin(); it != m_localPerturbations.end(); ++it) if (*it == i_perturbation) return true;
    return false;
}

void Spline::RemoveLocalPerturbation(LocalPerturbation *i_perturbationToRemove)
{
    for (std::vector<LocalPerturbation*>::iterator it = m_localPerturbations.begin(); it != m_localPerturbations.end(); ++it)
        if (*it == i_perturbationToRemove)
    {
        m_localPerturbations.erase(it);
        break;
    }
}

void Spline::getFitData(int &o_NData, Point *&o_dataPoints)
{
    o_NData = NDataPoints;
    o_dataPoints = DataPoints;
}

void Spline::GetSFuncs(int i, double** Data, int NPoints)
{
	bool DataEqual = false;
	int n;
	if (NPoints == NDataPoints)
	{
		for (n=0; n < NPoints && Data[n][0] == DataPoints[n].x; n++) ;
		if (n == NPoints) DataEqual = true;
	}
	if (!DataEqual)
	{
		delete[] DataPoints;
		DataPoints = new Point[NDataPoints = NPoints];
		for (n=0; n < NPoints; n++) DataPoints[n].x = Data[n][0];
		calcS();
	}
	for (n=0; n < NPoints; n++) Data[n][1] = S[i][n];
}

double Spline::gety(double x) const
{
	double DR, d6 = 1.0 / 6, DR2, A, B, C, D;
	int p=0;
	if (x < points[0].x || (NSplinePoints > 1 && points[0].x == points[1].x)) 
		return points[0].y;
	if (x > points[NSplinePoints - 1].x) return points[NSplinePoints - 1].y;
	while (x > points[p+1].x) p++;
	DR = points[p+1].x - points[p].x;
	DR2 = DR * DR * d6;
	A = (points[p+1].x - x) / DR;
	B = 1.0 - A;
	C = A * (A*A - 1.0) * DR2;
	D = B * (B*B - 1.0) * DR2;
	return A * points[p].y + B * points[p+1].y + C * points[p].yss + D * points[p+1].yss;
}

double Spline::getys(const double x) const
{
    int n;
    if (x < points[1].x) n=1;
    while (n < NSplinePoints - 1 && points[n].x < x) ++n;
    const double dx = points[n].x - points[n-1].x, ddx = 1.0 / dx, d6 = 1.0 / 6, A = (points[n].x - x) * ddx, B = 1.0 - A;
    return (points[n].y - points[n-1].y) * ddx - (3.0 * A*A - 1.0) * d6 * dx * points[n-1].yss + (3.0 * B*B - 1.0) * d6 * dx * points[n].yss;
}

double Spline::getyn(int n)
{
	if (Natural) return points[n].y;
	if (n==0) return points[n].yss;
	if (n == NSplinePoints + 1) return points[n-2].yss;
	return points[n-1].y;
}

bool Spline::readData(QTextStream *i_stream)
{
    bool endFound = false;
    int m;
    QStringList pointsList;
    while (!i_stream->atEnd())
    {
        QString Buffer = i_stream->readLine();
        if (Buffer.indexOf("IsNaturalSpline:", 0, Qt::CaseInsensitive) >= 0) Natural = (Buffer.indexOf("true", 0, Qt::CaseInsensitive) > 0);
        else if ((m = Buffer.indexOf("|")) > 0 && Buffer.left(m).toDouble() > 0.0) pointsList << Buffer;
        else if (Buffer.indexOf("End Spline", 0, Qt::CaseInsensitive) >= 0)
        {
            endFound = true;
            break;
        }
    }
    if (!endFound)
    {
        m_error = SEFileCorrupted;
        return false;
    }
    else if (pointsList.size() == 0)
    {
        m_error = SESplineEmpty;
        return false;
    }
    if (points != 0) delete[] points;
    points = new SplinePoint[pointsList.size()];
    for (int n=0; n < pointsList.size(); ++n)
    {
        QStringList currPoint = pointsList[n].split('|');
        if (currPoint.size() < 2)
        {
            m_error = SEFileCorrupted;
            return false;
        }
        points[n].x = currPoint[0].toDouble();
        points[n].y = currPoint[1].toDouble();
        points[n].yss = currPoint[2].toDouble();
    }
    NSplinePoints = pointsList.size();
    m_error = SENoError;
    return true;
}

void Spline::setNatural(bool natural)
{
	if (Natural != natural)
	{
		Natural = natural;
		deleteSL();
	}
}

void Spline::setSplinePoints(int NumSplinePoints, SplinePoint* SplinePoints)
{
	deleteSL();
	delete[] points;
	points = SplinePoints;
	NSplinePoints = NumSplinePoints;
}

void Spline::setyn(int n, double y)
{
	if (Natural) points[n].y = y;
	else if (n==0) points[n].yss = y;
	else if (n == NSplinePoints + 1) points[n-2].yss = y;
	else points[n-1].y = y;
}

void Spline::shifty(double val)
{
	int n;
	for (n=0; n < NSplinePoints; n++) points[n].y += val;
}

void Spline::UpdateY(double* C)
{
	int n, s = (Natural ? 0 : 1);
	if (!Natural) S0 += C[0];
	for (n=0; n < NSplinePoints; n++) points[n].y += C[n+s];
	if (!Natural) SN += C[n+s];
}

void Spline::writeData(QTextStream *i_stream) const
{
    *i_stream << "Begin Spline\n";
    *i_stream << "IsNaturalSpline: " << (Natural ? "true\n" : "false\n");
    *i_stream << "Column titles: x | y | d^2y/dx^2\n";
    for (int n=0; n < NSplinePoints; ++n) *i_stream << QString::number(points[n].x, 'f', 16) << " | " << QString::number(points[n].y, 'f', 16) << " | "
                                                    << QString::number(points[n].yss, 'f', 16) << "\n";
    *i_stream << "End Spline\n";
}

void Spline::MovePoint(const int i_pointIndex, const double i_newX, const double i_newY)
{
    points[i_pointIndex].x = i_newX;
    points[i_pointIndex].y = i_newY;
    if (S != 0)
    {
        Destroy(S, Natural ? NSplinePoints : NSplinePoints + 2);
        S = 0;
    }
    CalcLMatrix(L, points, 0.0, NSplinePoints, Natural);
    calcYss();
}

void Spline::RemovePoint(const int i_pointIndex)
{
    deleteSL();
    for (int i = i_pointIndex + 1; i < NSplinePoints; ++i) points[i-1] = points[i];
    NSplinePoints--;
    Refit();
}

void Spline::AddPoint(const double i_x, const double i_y)
{
    int i;
    deleteSL();
    SplinePoint* pointBuffer = new SplinePoint[NSplinePoints + 1];
    for (i=0; i < NSplinePoints && points[i].x < i_x; ++i) pointBuffer[i] = points[i];
    pointBuffer[i].x = i_x;
    pointBuffer[i].y = i_y;
    for (; i < NSplinePoints; ++i) pointBuffer[i+1] = points[i];
    delete[] points;
    points = pointBuffer;
    ++NSplinePoints;
    Refit();
}
