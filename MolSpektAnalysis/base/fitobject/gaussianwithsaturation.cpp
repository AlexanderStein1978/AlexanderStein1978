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


GaussianWithSaturation::GaussianWithSaturation(double *x, double *y, double *Sig, int N) : Gaussian(x, y, Sig, N)
{
}

GaussianWithSaturation::GaussianWithSaturation(const QString &data) : Gaussian(data)
{
}

GaussianWithSaturation::~GaussianWithSaturation()
{
}

double GaussianWithSaturation::GetPoint(double i_E) const
{
    return Offset - applySaturation(Offset - Gaussian::GetPoint(i_E));
}

void GaussianWithSaturation::getLineY(double *Ycalc) const
{
    Gaussian::getLineY(Ycalc);
    applySaturationOnLineProfile(Ycalc);
}

bool GaussianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    Gaussian::getCalcYAndDerivatives(Ycalc, deriv);
    applySaturationOnCalcYAndDerivatives(Ycalc, deriv);
    return true;
}

bool GaussianWithSaturation::getCalcY(double *Ycalc) const
{
    Gaussian::getLineY(Ycalc);
    applySaturationOnCalcY(Ycalc);
    return true;
}
