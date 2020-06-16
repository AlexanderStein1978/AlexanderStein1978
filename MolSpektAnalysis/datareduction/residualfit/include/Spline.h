//
// C++ Interface: Spline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SPLINE_H
#define SPLINE_H


#include "SplinePoint.h"

#include <vector>

class QTextStream;

class LocalPerturbation;

struct Point;


class Spline
{
public:

    enum SplineError{SENoError, SEFileCorrupted, SESplineEmpty};

    Spline();
	Spline(int NSplinePoints, SplinePoint* SPoints, bool natural);
    Spline(const Spline& i_other);
	virtual ~Spline();
    Spline& operator=(const Spline& i_right);
	void FitTo(int N, Point *points);
    void Refit();
    void JoinForDeperturbation(Spline* i_other, int JminPert, int JmaxPert);
    void RemoveLocalPerturbation(LocalPerturbation* i_perturbationToRemove);
    void MovePoint(const int i_pointIndex, const double i_newX, const double i_newY);
    void RemovePoint(const int i_pointIndex);
    void AddPoint(const double i_x, const double i_y);
    double gety(double x) const;
    double getys(const double x) const;
	double getyn(int n);
	void setyn(int n, double y);
	void calcYss();
	void GetSFuncs(int n, double **Data, int NPoints);
	void setNatural(bool natural);
	void setSplinePoints(int NumSplinePoints, SplinePoint* SplinePoints);
	void shifty(double val);
    bool readData(QTextStream *i_stream);
    void writeData(QTextStream *i_stream) const;
    void getFitData(int &o_NData, Point*& o_dataPoints);
    void addLocalPerturbation(LocalPerturbation* i_perturbation);
    bool containsLocalPerturbation(const LocalPerturbation* const i_perturbation) const;
    void copyData(const int N, const Point *const points);

    inline double getx0() const
	{
		return points[0].x;
	}
	
    inline double getxN() const
	{
		return points[NSplinePoints - 1].x;
	}
	
	inline int getNSFuncs()
	{
		return (NSplinePoints > 1 ? (Natural ? NSplinePoints : NSplinePoints + 2) : 0);
	}
	
	inline bool isNatural()
	{
		return Natural;
	}
	
	inline int getNPoints()
	{
		return NSplinePoints;
	}

    inline const SplinePoint& getSplinePoint(const int i) const
    {
        return points[i];
    }
	
	inline int getNPar()
	{
		return (Natural ? NSplinePoints : NSplinePoints + 2);
	}

    inline int getNumberOfLocalPerturbations() const
    {
        return m_localPerturbations.size();
    }

    inline const LocalPerturbation& getLocalPerturbation(int i_index) const
    {
        return *m_localPerturbations[i_index];
    }

    inline void RemoveLocalPerturbation(int i_index)
    {
        m_localPerturbations.erase(m_localPerturbations.begin() + i_index);
    }

    inline bool isFitDataAvailable() const
    {
        return (DataPoints != 0);
    }

    inline SplineError GetError()
    {
        return m_error;
    }
	
protected:
	void calcS();
	void calcFQSRes(double *Res, double &FQS);
	void UpdateY(double *C);
	void deleteSL();
    void copyFields(const Spline& i_other);
	
	double **L, **S, S0, SN;
	int NDataPoints, NSplinePoints;
    SplineError m_error;
	Point *DataPoints;
	SplinePoint *points;
    bool Natural;
    std::vector<LocalPerturbation*> m_localPerturbations;
};

#endif
