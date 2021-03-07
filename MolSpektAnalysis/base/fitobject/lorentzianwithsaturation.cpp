//
// C++ Implementation: Lorentzian with saturation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2021 - 2021
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

double LorentzianWithSaturation::GetPoint(double i_E) const
{
    return Offset - applySaturation(Offset - Lorentzian::GetPoint(i_E));
}


void LorentzianWithSaturation::getLineY(double *Ycalc) const
{
    Lorentzian::getLineY(Ycalc);
    applySaturationOnLineProfile(Ycalc);
}

bool LorentzianWithSaturation::getCalcYAndDerivatives(double *Ycalc, double **deriv)
{
    Lorentzian::getCalcY(Ycalc);
    Lorentzian::getDerivatives(deriv);
    applySaturationOnCalcYAndDerivatives(Ycalc, deriv);
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
    applySaturationOnCalcY(Ycalc);
    return true;
}

