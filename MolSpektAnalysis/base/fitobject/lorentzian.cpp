//
// C++ Implementation: Lorentzian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "lorentzian.h"
#include "gaussian.h"
#include "utils.h"
#include "tools.h"

#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>

#include <QFile>
#include <QTextStream>


using std::isnan;

Lorentzian::Lorentzian() : LineProfile(4), A(0.0), E(0.0), Gamma(0.0), AWF(0.0), EWF(0.0), GWF(0.0)
{
}

Lorentzian::Lorentzian(double* x, double* y, double* Sig, int N) : LineProfile(4)
{
	setData(x, y, Sig, N);
}

Lorentzian::Lorentzian(const QString &data) : LineProfile(4), A(0.0), E(0.0), Gamma(0.0), AWF(1.0), EWF(1.0), GWF(1.0)
{
    QStringList L = data.split('\t', QString::SkipEmptyParts);
    if (L.size() >= 7)
    {
        A = L[0].toDouble();
        E = L[1].toDouble();
        Gamma = L[2].toDouble();
        setOffset(L[3].toDouble());
        setDataRange(L[4].toDouble(), L[5].toDouble());
        setSubtracted(L[6] == "true");
    }
}

Lorentzian::Lorentzian(const Lorentzian &other) : LineProfile(other), A(other.A), E(other.E), Gamma(other.Gamma), AWF(other.AWF),
    EWF(other.EWF), GWF(other.GWF)
{
}

Lorentzian::Lorentzian(const Gaussian &other) : LineProfile(other), AWF(1.0), EWF(1.0), GWF(1.0)
{
    InitializeFromLineProfile(other);
}

Lorentzian::Lorentzian(const LineProfile &other) : LineProfile(other), AWF(1.0), EWF(1.0), GWF(1.0)
{
    InitializeFromLineProfile(other);
}

bool Lorentzian::getCalcY(double* YCalc) const
{
    int n;
	double AG = A * Gamma;
    for (n=0; n < nData; n++) YCalc[n] = Offset + AG / (sqr(E - X[n]) + Gamma);
    return true;
}

void Lorentzian::getLineY(double *Ycalc) const
{
    double AG = A * Gamma;
    for (int n=0; n < nData; n++) Ycalc[n] = AG / (sqr(E - X[n]) + Gamma);
}

double Lorentzian::GetPoint(double i_E) const
{
    return Offset + GetProfilePoint(E - i_E);
}

double Lorentzian::GetProfilePoint(double i_relE) const
{
    return A * Gamma / (sqr(i_relE) + Gamma);
}

void Lorentzian::getDerivatives(double** deriv)
{
    double ME = 0.0, MA = 0.0, MG = 0.0, D, DSq, Div, AG = A * Gamma;
	int n;
	for (n=0; n < nData; n++)
	{
		D = E - X[n];
		DSq = sqr(D);
		Div = 1.0 / (DSq + Gamma);
		deriv[n][0] = Gamma * Div;
        deriv[n][1] = -2.0 * AG * D * sqr(Div);
		deriv[n][2] = A * Div * (1.0 - Gamma * Div);
        deriv[n][3] = 1.0;
		if (fabs(deriv[n][0]) > MA) MA = fabs(deriv[n][0]);
		if (fabs(deriv[n][1]) > ME) ME = fabs(deriv[n][1]);
		if (fabs(deriv[n][2]) > MG) MG = fabs(deriv[n][2]);
	}
	AWF = 1.0 / MA;
	EWF = 1.0 / ME;
	GWF = 1.0 / MG;
	for (n=0; n < nData; n++)
	{
		deriv[n][0] *= AWF;
		deriv[n][1] *= EWF;
		deriv[n][2] *= GWF;
	}
}

void Lorentzian::getPar(double* Par) const
{
    Par[0] = A;
	Par[1] = E;
	Par[2] = Gamma;
    Par[3] = Offset;
}

void Lorentzian::GetValues(double &o_Intensity, double &o_CenterFreq, double &o_Width, double &o_Offset) const
{
    o_Intensity = A;
    o_CenterFreq = E;
    o_Width = 2.0 * sqrt(Gamma);
    o_Offset = Offset;
}

double Lorentzian::GetWidth() const
{
    return 2.0 * sqrt(Gamma);
}

void Lorentzian::InitializeFromLineProfile(const LineProfile &other)
{
    double width;
    other.GetValues(A, E, width, Offset);
    Gamma = sqr(0.5 * width);
}

void Lorentzian::Serialize(QTextStream &stream, const bool finish) const
{
    double Estart, Eend;
    GetDataRange(Estart, Eend);
    stream << "Lorentzian" << (isWithSaturation() ? "withSaturation\t" : "\t") << QString::number(A, 'g', 8) << '\t' << QString::number(E, 'g', 13) << '\t'
           << QString::number(Gamma, 'g', 8) << '\t' << QString::number(Offset, 'g', 8) << '\t' << QString::number(Estart, 'g', 13) << '\t' << QString::number(Eend, 'g', 13) << '\t'
           << (isLineSubtracted() ? "true" : "false");
    if (finish) stream << '\n';
}

void Lorentzian::setData(double* x, double* y, double* Sig, int N)
{
    int n, m=0;
	AWF = EWF = GWF = 0.0;
	for (n=0, A=0.0; n<N; n++) if (y[n] > A)
	{
		m=n;
		A = y[n];
	}
	E = x[m];
	if (N>2 && m>0 && m<N-1) Gamma = sqr(x[m+1] - x[m-1]);
	else Gamma = 1;
    if (N > 0)
    {
        Offset = x[0];
        setDataRange(x[0], x[N-1]);
    }
	FitObject::setData(x, y, Sig, N);
}

void Lorentzian::setPar(double* Par)
{
    A = Par[0];
	E = Par[1];
	Gamma = Par[2];
    Offset = Par[3];
}

void Lorentzian::SetValues(const double Intensity, const double CenterFreq, const double Width, const double iOffset)
{
    A = Intensity;
    E = CenterFreq;
    Gamma = sqr(0.5 * Width);
    Offset = iOffset;
}

void Lorentzian::updatePar(double* C)
{
    A += AWF * C[0];
	E += EWF * C[1];
	Gamma += GWF * C[2];
    Offset += C[3];
}
