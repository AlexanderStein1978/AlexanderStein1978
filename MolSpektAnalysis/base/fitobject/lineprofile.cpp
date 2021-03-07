//
// C++ Implementation: LineProfile
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2021 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "lineprofile.h"


LineProfile::LineProfile() : Offset(0.0), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
}

LineProfile::LineProfile(int nPar) : FitObject(nPar), Offset(0.0), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
}

LineProfile::LineProfile(int NPar, double *x, double *y, double *Sig, int N) : FitObject(NPar, x, y, Sig, N), Offset(0.0), isSubtracted(false), m_Estart(0.0), m_Eend(0.0)
{
    if (N>0)
    {
        m_Estart = x[0];
        m_Eend = x[N-1];
    }
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

void LineProfile::setDataRange(const double E_start, const double E_end)
{
    m_Estart = E_start;
    m_Eend = E_end;
}
