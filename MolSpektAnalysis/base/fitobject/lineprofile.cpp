//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "lineprofile.h"


LineProfile::LineProfile() : Offset(0.0), hideSaturation(false), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
}

LineProfile::LineProfile(int nPar) : FitObject(nPar), Offset(0.0), hideSaturation(false), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
}

LineProfile::LineProfile(int NPar, double *x, double *y, double *Sig, int N) : FitObject(NPar, x, y, Sig, N), Offset(0.0),
    hideSaturation(false), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
    if (N>0)
    {
        m_Estart = x[0];
        m_Eend = x[N-1];
    }
}

LineProfile::LineProfile(const LineProfile &other) : FitObject(other), Offset(other.Offset), hideSaturation(false),
    isSubtracted(other.isSubtracted), m_Estart(other.m_Estart), m_Eend(other.m_Eend)
{
}

void LineProfile::applySaturationOnLineProfile(double *Ycalc) const
{
    for (int n=0; n < nData; ++n) Ycalc[n] = -applySaturation(-Ycalc[n]);
}

void LineProfile::applySaturationOnCalcYAndDerivatives(double *Ycalc, double **deriv) const
{
    double OSq = Offset * Offset;
    for (int n=0; n < nData; n++)
    {
        double I(Offset - Ycalc[n]), sum(Offset + I), sf = OSq / (sum * sum);
        Ycalc[n] = Offset - applySaturation(I);
        for (int i=0; i < nPar; ++i) deriv[n][i] *= sf;
    }
}

void LineProfile::applySaturationOnCalcY(double *Ycalc) const
{
    for (int n=0; n < nData; ++n) Ycalc[n] = Offset - applySaturation(Offset - Ycalc[n]);
}

void LineProfile::GetDataRange(double &o_Estart, double &o_Eend) const
{
    o_Estart = m_Estart;
    o_Eend = m_Eend;
}

void LineProfile::GetDataRange(double &Emin, double &Imin, double &Emax, double &Imax) const
{
    GetDataRange(Emin, Emax);
    Imin = 1e99;
    Imax = -1e99;
    for (int n=0; n < nData; ++n)
    {
        if (Y[n] < Imin) Imin = Y[n];
        else if (Y[n] > Imax) Imax = Y[n];
    }
}

QString LineProfile::getProfileTypeName(const LineProfileType type)
{
    if (type == GaussianType) return "Gaussian";
    return "Lorentzian";
}

void LineProfile::setDataRange(const double E_start, const double E_end)
{
    m_Estart = E_start;
    m_Eend = E_end;
}
