//
// C++ Implementation: FitObject
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2016
//
// Copyright: See README file that comes with this source code
//
//


#include "fitobject.h"
#include "utils.h"
#include "tools.h"

#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>

#include <QFile>
#include <QTextStream>


using std::isnan;


FitObject::FitObject() : X(0), Y(0), sig(0), nData(0), nPar(0), m_maxBadBest(1), m_derivativesNotCalculated(false), m_DebugLogFile(0), m_DebugLogStream(0)
{
}

FitObject::FitObject(int NPar) : m_maxBadBest(1), m_derivativesNotCalculated(false), m_DebugLogFile(0), m_DebugLogStream(0)
{
	nData = 0;
	nPar = NPar;
	X = Y = sig = 0;
}

FitObject::FitObject(int NPar, double* x, double* y, double* Sig, int N) : m_maxBadBest(1), m_derivativesNotCalculated(false), m_DebugLogFile(0),
    m_DebugLogStream(0)
{
	nPar = NPar;
	X=x;
	Y=y;
	sig = Sig;
	nData = N;
}

FitObject::FitObject(const FitObject &i_toCopy) : m_maxBadBest(1)
{
    copyData(i_toCopy);
}

FitObject::~FitObject()
{
	destroyData();
}

FitObject& FitObject::operator=(const FitObject& i_right)
{
    if (this != &i_right)
    {
        destroyData();
        copyData(i_right);
    }
    return *this;
}

void FitObject::addDataPoint(double x, double y, double Sig)
{
	double *XB = new double[++nData], *YB = new double[nData], *SB = new double[nData];
	int n;
	for (n=0; (n < nData - 1 ? X[n] < x : false); n++)
	{
		XB[n] = X[n];
		YB[n] = Y[n];
		SB[n] = sig[n];
	}
	XB[n] = x;
	YB[n] = y;
	SB[n] = Sig;
	while (n+1 < nData)
	{
		XB[n+1] = X[n];
		YB[n+1] = Y[n];
		SB[n+1] = sig[n];
		n++;
	}
	destroyData();
	X = XB;
	Y = YB;
	sig = SB;
	
}

void FitObject::copyData(const FitObject &i_objectToCopyDataFrom)
{
    nPar = i_objectToCopyDataFrom.nPar;
    nData = i_objectToCopyDataFrom.nData;
    X = new double[nData];
    Y = new double[nData];
    sig = new double[nData];
    m_derivativesNotCalculated = i_objectToCopyDataFrom.m_derivativesNotCalculated;
    m_DebugLogFile = i_objectToCopyDataFrom.m_DebugLogFile;
    m_DebugLogStream = i_objectToCopyDataFrom.m_DebugLogStream;
    for (int n=0; n < nData; ++n)
    {
        X[n] = i_objectToCopyDataFrom.X[n];
        Y[n] = i_objectToCopyDataFrom.Y[n];
        sig[n] = i_objectToCopyDataFrom.sig[n];
    }

}

void FitObject::destroyData()
{
	if (X != 0)
	{
		delete[] X;
		delete[] Y;
		delete[] sig;
	}
}

bool FitObject::getCalcYAndDerivatives(double *Ycalc, double **/*deriv*/)
{
    m_derivativesNotCalculated = true;
    return getCalcY(Ycalc);
}

bool FitObject::getCalcY(double* Ycalc) const
{
	int n;
	for (n=0; n < nData; n++) Ycalc[n] = 0.0;
    return false;
}

void FitObject::getDerivatives(double** deriv)
{
	int n, m;
	for (n=0; n < nData; n++) for (m=0; m < nPar; m++) deriv[n][m] = 0.0;
}

void FitObject::getPar(double* Par)
{
	int n;
	for (n=0; n < nPar; n++) Par[n] = 0;
}

void FitObject::setNPar()
{
}

double FitObject::GetSigma() const
{
    if (nData <= 0) return -1.0;
    double Sigma, FQS = 0.0, *YC = new double[nData];
    bool result = getCalcY(YC);
    if (result)
    {
        for (int n=0; n < nData; ++n) FQS += sqr(Y[n] - YC[n]) * sig[n];
        Sigma = sqrt(FQS / (nData - 1.0));
    }
    else Sigma = -1.0;
    delete[] YC;
    return Sigma;
}

double FitObject::LevenbergMarquardt(int MaxIt, double MinImp)
{
    setNPar();
    if (nPar <= 0) return -1.0;
	double **LMEQS = Create(nPar, nPar + 1), **tLMEQS = Create(nPar, nPar + 1), wt, **EQS = Create(nData, nPar);
	double *YC = new double[nData], done, BadCount, lambda, *C = new double[nPar], FQS, aFQS, FQSv, *bPar = new double[nPar];
    bool badBest = false;
    int n, i, j, it = 0, badBestCounter = 0;
	for (n=0; n < nPar; n++) for (i=0; i <= nPar; i++) LMEQS[n][i] = 0.0;
    getCalcYAndDerivatives(YC, EQS);
    if (m_derivativesNotCalculated) getDerivatives(EQS);
	getPar(bPar);
	for (FQSv = 0.0, n=0; n < nData; n++) FQSv += sqr(Y[n] - YC[n]) * sig[n];
	for (n=0; n < nData; n++)
	{
		for (i=0; i < nPar; i++)
		{
			wt = EQS[n][i] * sig[n];
			for (j=0; j < nPar; j++) LMEQS[i][j] += wt * EQS[n][j];
			LMEQS[i][nPar] += wt * (Y[n] - YC[n]);
		}
	}
	aFQS = FQS = FQSv;
    for (done = BadCount = 0, lambda = 1e-3; true; it++)
	{
		//printf("bPar[0]=%g, bPar[1]=%g, bPar[2]=%g, aFQS=%g, FQS=%g\n", bPar[0], bPar[1], bPar[2], aFQS, FQS);
		if (done == 4) lambda = 0.0;
		for (n=0; n < nPar; n++)
		{
			for (i=0; i <= nPar; i++) tLMEQS[n][i] = LMEQS[n][i];
			tLMEQS[n][n] *= (1.0 + lambda);
		}
		SolvLinEqS(tLMEQS, nPar);
		for (n=0; n < nPar; n++) C[n] = tLMEQS[n][nPar];
		updatePar(C);
        bool goodParSet = getCalcYAndDerivatives(YC, EQS);
		for (FQS = 0.0, n=0; n < nData; n++) FQS += sqr(Y[n] - YC[n]) * sig[n];
        if (0 != m_DebugLogStream)
            *m_DebugLogStream << QString("FQS=%1\tFQSv=%2\taFQS=%3\n").arg(FQS, 0, 'g', 6).arg(FQSv, 0, 'g', 6).arg(aFQS, 0, 'g', 6);
        if (!goodParSet && FQS < aFQS)
        {
            badBest = true;
            ++badBestCounter;
        }
        else badBest = false;
        if (it == MaxIt || done == 4 || BadCount == 20 || badBestCounter == m_maxBadBest) break;
		if (abs(FQS - aFQS) <  1e-10 * aFQS || (FQS < aFQS && aFQS - FQS < MinImp)) done++;
        if (FQS < aFQS && !isnan(FQS) && !badBest)
		{
			BadCount = 0;
            badBestCounter = 0;
			lambda *= 0.1;
			aFQS = FQS;
            if (m_derivativesNotCalculated) getDerivatives(EQS);
			getPar(bPar);
			for (n=0; n < nPar; n++) for (i=0; i <= nPar; i++) LMEQS[n][i] = 0.0;
			for (n=0; n < nData; n++) for (i=0; i < nPar; i++)
			{
				wt = EQS[n][i] * sig[n];
				for (j=0; j < nPar; j++) LMEQS[i][j] += wt * EQS[n][j];
				LMEQS[i][nPar] += wt * (Y[n] - YC[n]);
			}
		}
		else 
		{
			BadCount++;
			lambda *= 10.0;
			setPar(bPar);
		}
	}
    if (FQS > aFQS || isnan(FQS) || badBest)
	{
		setPar(bPar);
		FQS = aFQS;
		//getCalcY(YC);
	}
	//for (n=0; n < nData; n++) printf("X[%d]=%g, YC=%g, YO=%g\n", n, X[n], YC[n], Y[n]);
	Destroy(LMEQS, nPar);
	Destroy(tLMEQS, nPar);
	Destroy(EQS, nData);
	delete[] YC;
	delete[] C;
	delete[] bPar;
	return FQS;
}

void FitObject::setData(double* x, double* y, double* Sig, int N)
{
	destroyData();
	X = x;
	Y = y;
	nData = N;
	sig = Sig;
}

void FitObject::setPar(double* /*Par*/)
{

}

void FitObject::updatePar(double* /*C*/)
{

}

void FitObject::InitDebugLogging(QString i_FileName)
{
    m_DebugLogFile = new QFile(i_FileName);
    m_DebugLogFile->open(QIODevice::WriteOnly);
    m_DebugLogStream = new QTextStream(m_DebugLogFile);
}

void FitObject::EndDebugLogging()
{
    delete m_DebugLogStream;
    m_DebugLogStream = 0;
    delete m_DebugLogFile;
    m_DebugLogFile = 0;
}
