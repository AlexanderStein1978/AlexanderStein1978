//
// C++ Implementation: Gaussian line profile with saturation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "gaussianwithsaturation.h"

#include <QStringList>
#include <QTextStream>

#include <cmath>


GaussianWithSaturation::GaussianWithSaturation(double *x, double *y, double *Sig, int N) : Gaussian(x, y, Sig, N), SB(B)
{
    double ym(y[0]);
    int nm;
    for (int n=0; n<N; ++n) if (y[n] < ym)
    {
        ym = y[n];
        nm = n;
    }
    IBo = 0.5 * (y[0] + y[nm]);
    iBn = getNBorderByI(IBo, nm, -1);
    oBn = getNBorderByI(IBo, nm, 1);
}

GaussianWithSaturation::GaussianWithSaturation(const QString &data) : Gaussian(), SB(0.0), IBo(0.0), iBn(0), oBn(0)
{
    QStringList L = data.split('\t', QString::SkipEmptyParts);
    if (L.size() >= 9)
    {
        initialize(L);
        SB = L[7].toDouble();
        IBo = L[8].toDouble();
    }
}

GaussianWithSaturation::~GaussianWithSaturation()
{
}

void GaussianWithSaturation::GetValues(double &o_Intensity, double &o_SatIntensity, double &o_SatBorder, double &o_CenterFreq, double &o_Width, double &o_Offset) const
{
    Gaussian::GetValues(o_Intensity, o_CenterFreq, o_Width, o_Offset);
    o_SatIntensity = SB;
    o_SatBorder = IBo;
}

void GaussianWithSaturation::SetValues(const double Intensity, const double SatIntensity, const double SatBorder, const double CenterFreq, const double Width, const double Offset)
{
    Gaussian::SetValues(Intensity, CenterFreq, Width, Offset);
    SB = SatIntensity;
    IBo = SatBorder;
}

double GaussianWithSaturation::GetPoint(double i_E) const
{
    if (nullptr != X)
    {
        if (i_E < X[iBn] || i_E > X[oBn]) return Gaussian::GetPoint(i_E);
        return GetPoint(i_E, LineSaturated);
    }
}

double GaussianWithSaturation::GetProfilePoint(double i_E) const
{
    if (nullptr != X)
    {
        if (i_E < X[iBn] || i_E > X[oBn]) return Gaussian::GetProfilePoint(i_E);
        double arg = (i_E - E) / G;
        return Offset * (1.0 - SB * exp(arg * arg));
    }
}

void GaussianWithSaturation::Serialize(QTextStream &stream, const bool finish) const
{
    Gaussian::Serialize(stream, false);
    stream << '\t' << QString::number(SB, 'g', 8) << '\t' << QString::number(IBo, 'g', 8);
    if (finish) stream << '\n';
}

void GaussianWithSaturation::getLineY(double *Ycalc) const
{
    double dG = 1.0 / G;
    for (int n=0; n < nData; ++n)
    {
        double arg = (X[n] - E) * dG;
        Ycalc[n] = (n < iBn || n > oBn ? B * exp(-arg * arg) : Offset * (1.0 - SB * exp(arg * arg)));
    }
}

void GaussianWithSaturation::updateBorders()
{
    iBn = findBorder(iBn, -1);
    oBn = findBorder(oBn, 1);
    IBo = (abs(GetPoint(X[iBn], LineGaussian) - GetPoint(X[iBn], LineSaturated)) < abs(GetPoint(X[oBn], LineGaussian) - GetPoint(X[oBn], LineSaturated)) ? Y[iBn] : Y[oBn]);
}

void GaussianWithSaturation::setData(double *x, double *y, double *Sig, int N)
{
    setData(x, y, Sig, N);
    iBn = getNBorderByI(IBo, (iBn > 0 && iBn < N - 1 ? iBn : 1), -1);
    oBn = getNBorderByI(IBo, (oBn > 0 && oBn < N - 1 ? oBn : iBn), 1);
}

bool GaussianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    int n;
    double dG = 1.0 / G;
    for (n=0; n < nData; n++)
    {
        bool satPart = n > iBn && n < oBn;
        double arg = (X[n] - E) * dG;
        double ex = exp((satPart ? arg : -arg) * arg);
        Ycalc[n] = (satPart ? Offset * SB * ex : Offset + B * ex);
        deriv[n][0] = (satPart ? 0.0 : ex);
        deriv[n][1] = dG * 2.0 * (satPart ? SB : B) * ex * arg;;
        if (satPart) deriv[n][1] *= (-Offset);
        deriv[n][2] = arg * deriv[n][1];
        deriv[n][3] = (satPart ? SB * ex : 1.0);
        deriv[n][4] = (satPart ? Offset * ex : 0.0);
    }
    return true;
}

bool GaussianWithSaturation::getCalcY(double *Ycalc) const
{
    double dG = 1.0 / G;
    for (int n=0; n < nData; ++n)
    {
        double arg = (X[n] - E) * dG;
        Ycalc[n] = (n <= iBn || n >= oBn ? Offset + B * exp(-arg * arg) : Offset * SB * exp(arg * arg));
    }
    return true;
}

void GaussianWithSaturation::getPar(double *Par)
{
    Gaussian::getPar(Par);
    Par[4] = SB;
}

void GaussianWithSaturation::setPar(double *Par)
{
    Gaussian::setPar(Par);
    SB = Par[4];
    updateBorders();
}

void GaussianWithSaturation::updatePar(double *C)
{
    Gaussian::updatePar(C);
    SB += C[4];
    updateBorders();
}

int GaussianWithSaturation::getNBorderByI(const double bY, const int nStart, const int sDirection) const
{
    int n = nStart;
    while (n > 0 && n < nData - 1 && Y[n] < bY) n += sDirection;
    while (n > 0 && n < nData - 1 && Y[n] > bY) n -= sDirection;
    if (n + sDirection > 0 && n + sDirection < nData - 1 && abs(Y[n + sDirection] - bY) < abs(Y[n] - bY)) n += sDirection;
    return n;
}

int GaussianWithSaturation::findBorder(const int nStart, const int sDirection) const
{
    int n = nStart;
    while (n > 0 && n < nData - 1 && GetPoint(X[n], LineGaussian) < GetPoint(X[n], LineSaturated)) n += sDirection;
    while (n > 0 && n < nData - 1 && GetPoint(X[n], LineGaussian) > GetPoint(X[n], LineSaturated)) n -= sDirection;
    if (n + sDirection > 0 && n + sDirection < nData - 1
        && abs(GetPoint(X[n + sDirection], LineGaussian) - GetPoint(X[n + sDirection], LineSaturated))
            < abs(GetPoint(X[n], LineGaussian) - GetPoint(X[n], LineSaturated))) n += sDirection;
    return n;
}

double GaussianWithSaturation::GetPoint(double i_E, const WhichLine l) const
{
    switch(l)
    {
    case LineBestfit:
        return GetPoint(i_E);
    case LineGaussian:
        return Gaussian::GetPoint(i_E);
    case LineSaturated:
        double arg = (i_E - E) / G;
        return Offset * SB * exp(arg * arg);
        break;
    }
}
