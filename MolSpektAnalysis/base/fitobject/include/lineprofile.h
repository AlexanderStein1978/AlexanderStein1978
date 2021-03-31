//
// C++ Interface: Line profile
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2021 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LINEPOFILE_H
#define LINEPOFILE_H


#include "fitobject.h"

class QString;


class LineProfile : public FitObject
{
public:
    enum LineProfileType
    {
        GaussianType,
        LorentzianType
    };

    LineProfile();
    LineProfile(int nPar);
    LineProfile(int NPar, double *x, double *y, double *Sig, int N);
    LineProfile(const LineProfile& other);

    void GetDataRange(double& o_Estart, double& o_Eend) const;
    void GetDataRange(double& Emin, double&Imin, double &Emax, double& Imax) const;
    virtual double GetPoint(double i_E) const = 0;
    virtual double GetProfilePoint(double i_relE) const = 0;
    virtual void Serialize(QTextStream& stream, const bool finish = true) const = 0;
    virtual void getLineY(double *Ycalc) const = 0;
    virtual double GetECenter() const = 0;
    virtual bool isValid() const = 0;
    virtual bool isWithSaturation() const = 0;
    virtual LineProfileType getType() const = 0;
    virtual void GetValues(double& o_Intensity, double& o_CenterFreq, double& o_Width, double& o_Offset) const = 0;
    virtual void SetValues(const double Intensity, const double CenterFreq, const double Width, const double Offset) = 0;

    inline void setSubtracted(const bool newValue)
    {
        isSubtracted = newValue;
    }

    inline bool isLineSubtracted() const
    {
        return isSubtracted;
    }

    inline void setOffset(const double newOffset)
    {
        Offset = newOffset;
    }

    inline double getOffset() const
    {
        return Offset;
    }

    inline void setHideSaturation(const bool hide)
    {
        hideSaturation = hide;
    }

protected:
    double Offset;
    bool hideSaturation, isSubtracted;

    void setDataRange(const double E_start, const double E_end);
    void applySaturationOnLineProfile(double *Ycalc) const;
    void applySaturationOnCalcYAndDerivatives(double *Ycalc, double **deriv) const;
    void applySaturationOnCalcY(double *Ycalc) const;

    inline double applySaturation(const double I) const
    {
        return Offset * I / (Offset + I);
    }

private:
    double m_Estart, m_Eend;
};

#endif
