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
    for (int n=0; n < nData; ++n) Ycalc[n] = -applySaturation(-Ycalc[n]);
}

bool GaussianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    Gaussian::getCalcYAndDerivatives(Ycalc, deriv);
    double OSq = Offset * Offset;
    for (int n=0; n < nData; n++)
    {
        double I(Offset - Ycalc[n]), sum(Offset + I), sf = OSq / (sum * sum);
        Ycalc[n] = Offset - applySaturation(I);
        for (int i=0; i<4; ++i) deriv[n][i] *= sf;
    }
    return true;
}

bool GaussianWithSaturation::getCalcY(double *Ycalc) const
{
    Gaussian::getLineY(Ycalc);
    for (int n=0; n < nData; ++n) Ycalc[n] = Offset - applySaturation(-Ycalc[n]);
    return true;
}
