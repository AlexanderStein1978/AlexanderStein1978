//
// C++ Interface: FitObject
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FITOBJECT_H
#define FITOBJECT_H


#include <QString>


class QFile;
class QTextStream;


class FitObject
{
public:
    FitObject();
	FitObject(int NPar);
	FitObject(int NPar, double *x, double *y, double *Sig, int N);
    FitObject(const FitObject& i_toCopy);
	~FitObject();
    FitObject& operator=(const FitObject& i_right);
	void addDataPoint(double x, double y, double Sig);
	virtual double LevenbergMarquardt(int MaxIt, double MinImp);
	virtual void setData(double *x, double *y, double *Sig, int N);
    void InitDebugLogging(QString FileName);
    void EndDebugLogging();
    double GetSigma() const;

    inline int GetNData() const
    {
        return nData;
    }
		
protected:
    virtual bool getCalcYAndDerivatives(double *Ycalc, double **deriv);
	virtual void getDerivatives(double **deriv);
    virtual bool getCalcY(double *Ycalc) const;
    virtual void getPar(double *Par) const = 0;
    virtual void setPar(double *Par) = 0;
    virtual void updatePar(double *C) = 0;
    virtual void setNPar() = 0;
	
	double *X, *Y, *sig;
    int nData, nPar, m_maxBadBest;

    bool m_derivativesNotCalculated;

    QFile* m_DebugLogFile;
    QTextStream* m_DebugLogStream;
	
private:
	void destroyData();
    void copyData(const FitObject& i_objectToCopyDataFrom);
};

#endif // FITOBJECT_H
