//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "lorentzianwithsaturation.h"


LorentzianWithSaturation::LorentzianWithSaturation()
{
}

LorentzianWithSaturation::LorentzianWithSaturation(double *x, double *y, double *Sig, int N) : Lorentzian(x, y, Sig, N)
{
}

LorentzianWithSaturation::LorentzianWithSaturation(const QString& data) : Lorentzian(data)
{
}

LorentzianWithSaturation::LorentzianWithSaturation(const Lorentzian &other) : Lorentzian(other)
{
}

LorentzianWithSaturation::LorentzianWithSaturation(const Gaussian &other) : Lorentzian(other)
{
}

LorentzianWithSaturation::LorentzianWithSaturation(const LineProfile &other) : Lorentzian(other)
{
}

double LorentzianWithSaturation::GetPoint(double i_E) const
{
    if (hideSaturation) return Lorentzian::GetPoint(i_E);
    return Offset - applySaturation(Offset - Lorentzian::GetPoint(i_E));
}


void LorentzianWithSaturation::getLineY(double *Ycalc) const
{
    Lorentzian::getLineY(Ycalc);
    if (hideSaturation) return;
    applySaturationOnLineProfile(Ycalc);
}

bool LorentzianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    Lorentzian::getCalcY(Ycalc);
    Lorentzian::getDerivatives(deriv);
    if (!hideSaturation) applySaturationOnCalcYAndDerivatives(Ycalc, deriv);
    return true;
}

void LorentzianWithSaturation::getDerivatives(double **deriv)
{
    double *Ycalc = new double[nData];
    getCalcYAndDerivatives(Ycalc, deriv);
    delete[] Ycalc;
}

bool LorentzianWithSaturation::getCalcY(double *Ycalc) const
{
    Lorentzian::getCalcY(Ycalc);
    if (!hideSaturation) applySaturationOnCalcY(Ycalc);
    return true;
}
