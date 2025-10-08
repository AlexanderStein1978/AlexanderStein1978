//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "gaussianwithsaturation.h"

#include <QStringList>
#include <QTextStream>

#include <cmath>


GaussianWithSaturation::GaussianWithSaturation(double *x, double *y, double *Sig, int N) : Gaussian(x, y, Sig, N)
{
}

GaussianWithSaturation::GaussianWithSaturation(const QString &data) : Gaussian(data)
{
}

GaussianWithSaturation::GaussianWithSaturation(const Gaussian& other) : Gaussian(other)
{
}

GaussianWithSaturation::GaussianWithSaturation(const Lorentzian& other) : Gaussian(other)
{
}

GaussianWithSaturation::GaussianWithSaturation(const LineProfile& other) : Gaussian(other)
{
}

GaussianWithSaturation::~GaussianWithSaturation()
{
}

double GaussianWithSaturation::GetPoint(double i_E) const
{
    if (hideSaturation) return Gaussian::GetPoint(i_E);
    return Offset - applySaturation(Offset - Gaussian::GetPoint(i_E));
}

void GaussianWithSaturation::getLineY(double *Ycalc) const
{
    Gaussian::getLineY(Ycalc);
    if (hideSaturation) return;
    applySaturationOnLineProfile(Ycalc);
}

bool GaussianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    Gaussian::getCalcYAndDerivatives(Ycalc, deriv);
    if (!hideSaturation) applySaturationOnCalcYAndDerivatives(Ycalc, deriv);
    return true;
}

bool GaussianWithSaturation::getCalcY(double *Ycalc) const
{
    Gaussian::getLineY(Ycalc);
    if (!hideSaturation) applySaturationOnCalcY(Ycalc);
    return true;
}
