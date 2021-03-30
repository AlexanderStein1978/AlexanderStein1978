//
// C++ Implementation: Gaussian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "gaussian.h"
#include "lorentzian.h"
#include "utils.h"
#include "tools.h"

#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>

#include <QFile>
#include <QTextStream>


using std::isnan;


Gaussian::Gaussian() : LineProfile(4), B(0.0), E(0.0), G(0.0)
{
}

Gaussian::Gaussian(double *x, double *y, double *Sig, int N) : LineProfile(4, x, y, Sig, N)
{
    Offset = y[0];
    double min = Offset, max = Offset, mx = 0.0, Mx = 0.0;
    int mn = 0, Mn = 0;
    for (int n=0; n < N; ++n)
    {
        if (y[n] < min)
        {
            min = y[n];
            mn = n;
            mx = x[n];
        }
        else if (y[n] > max)
        {
            max = y[n];
            Mn = n;
            Mx = x[n];
        }
    }
    double HI, FWHM;
    int cn;
    if (max - Offset > Offset - min)
    {
        E = Mx;
        B = max - Offset;
        HI = 0.5 * (max + Offset);
        cn = Mn;
    }
    else
    {
        E = mx;
        B = min - Offset;
        HI = 0.5 * (min + Offset);
        cn = mn;
    }
    int n;
    for (n = cn; n<N && (B > 0.0 ? y[n] > HI : y[n] < HI); ++n) ;
    FWHM = (abs(y[n] - HI) < abs(y[n-1] - HI) ? x[n] : x[n-1]);
    for (n = cn; n>0 && (B > 0.0 ? y[n] > HI : y[n] < HI); --n) ;
    FWHM -= (abs(y[n] - HI) < abs(y[n+1] - HI) ? x[n] : x[n+1]);
    G = 0.6005 * FWHM;
}

Gaussian::Gaussian(const QString &data) : LineProfile(4), B(0.0), E(0.0), G(0.0)
{
    QStringList L = data.split('\t', QString::SkipEmptyParts);
    initialize(L);
}

Gaussian::Gaussian(const Gaussian &other) : LineProfile(other), B(other.B), E(other.E), G(other.G)
{
}

Gaussian::Gaussian(const Lorentzian &other) : LineProfile(other), B(other.GetAmplitude()), E(other.GetECenter()),
    G(0.6005 * other.GetWidth())
{
}

Gaussian::Gaussian(const LineProfile &other) : LineProfile(other)
{
    double width;
    other.GetValues(B, E, width, Offset);
    G = 0.6005 * width;
}

Gaussian::~Gaussian()
{
}

void Gaussian::initialize(const QStringList &data)
{
    if (data.size() >= 7)
    {
        B = data[0].toDouble();
        E = data[1].toDouble();
        G = data[2].toDouble();
        setOffset(data[3].toDouble());
        setDataRange(data[4].toDouble(), data[5].toDouble());
        setSubtracted(data[6] == "true");
    }
}

void Gaussian::GetValues(double &o_B, double &o_E, double &o_Width, double &o_Offset) const
{
    o_B = B;
    o_E = E;
    o_Width = 1.665 * G;
    o_Offset = getOffset();
}

void Gaussian::SetValues(const double Intensity, const double CenterFreq, const double Width, const double iOffset)
{
    B = Intensity;
    E = CenterFreq;
    G = 0.6006 * Width;
    Offset = iOffset;
}

bool Gaussian::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    int n;
    double dG = 1.0 / G;
    for (n=0; n < nData; n++)
    {
        double arg = (X[n] - E) * dG;
        double ex = exp(-arg * arg);
        double fac = 2.0 * B * ex * arg;
        Ycalc[n] = Offset + B * ex;
        deriv[n][0] = ex;
        deriv[n][1] = dG * fac;
        deriv[n][2] = arg * deriv[n][1];
        deriv[n][3] = 1.0;
    }
    return true;
}

bool Gaussian::getCalcY(double *Ycalc) const
{
    double dG = 1.0 / G;
    for (int n=0; n < nData; ++n)
    {
        double arg = (X[n] - E) * dG;
        Ycalc[n] = Offset + B * exp(-arg * arg);
    }
    return true;
}

void Gaussian::getLineY(double *Ycalc) const
{
    double dG = 1.0 / G;
    for (int n=0; n < nData; ++n)
    {
        double arg = (X[n] - E) * dG;
        Ycalc[n] = B * exp(-arg * arg);
    }
}

double Gaussian::GetPoint(double i_E) const
{
    double arg = (i_E - E) / G;
    return Offset + B * exp(-arg * arg);
}

double Gaussian::GetProfilePoint(double i_relE) const
{
    double arg = i_relE / G;
    return exp(-arg * arg);
}

void Gaussian::getPar(double *Par) const
{
    Par[0] = B;
    Par[1] = E;
    Par[2] = G;
    Par[3] = Offset;
}

void Gaussian::setPar(double *Par)
{
    B = Par[0];
    E = Par[1];
    G = Par[2];
    Offset = Par[3];
}

void Gaussian::updatePar(double *C)
{
    B += C[0];
    E += C[1];
    G += C[2];
    Offset += C[3];
}

void Gaussian::Serialize(QTextStream &stream, const bool finish) const
{
    double Estart, Eend;
    GetDataRange(Estart, Eend);
    stream << "Gaussian" << (isWithSaturation() ? "withSaturation\t" : "\t") << QString::number(B, 'g', 8) << '\t' << QString::number(E, 'g', 13) << '\t' << QString::number(G, 'g', 8)
           << '\t' << QString::number(Offset, 'g', 8) << "\t" << QString::number(Estart, 'g', 13) << '\t' << QString::number(Eend, 'g', 13) << '\t'
           << (isLineSubtracted() ? "true" : "false");
    if (finish) stream << '\n';
}
