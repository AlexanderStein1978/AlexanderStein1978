//
// C++ Interface: Lorentzian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LORENTZIAN_H
#define LORENTZIAN_H


#include "fitobject.h"

class QString;


class Gaussian : public FitObject
{
public:
    Gaussian(double *x, double *y, double *Sig, int N);
    Gaussian(const QString& data);
    virtual ~Gaussian();
    void GetValues(double& o_Intensity, double& o_CenterFreq, double& o_Width, double& o_Offset) const;
    void SetValues(const double Intensity, const double CenterFreq, const double Width, const double Offset);
    void GetDataRange(double& o_Estart, double& o_Eend) const;
    double GetPoint(double i_E) const;
    double GetProfilePoint(double i_relE) const;
    void Serialize(QTextStream& stream) const;
    void getLineY(double *Ycalc) const;

    inline double GetEcenter() const
    {
        return E;
    }

    inline bool isValid() const
    {
        return B != 0.0;
    }

    inline void setSubtracted(const bool newValue)
    {
        isSubtracted = newValue;
    }

    inline bool isLineSubtracted() const
    {
        return isSubtracted;
    }

protected:
    bool getCalcYAndDerivatives(double *Ycalc, double **deriv) override;
    bool getCalcY(double *Ycalc) const override;
    void getPar(double *Par) override;
    void setPar(double *Par) override;
    void updatePar(double *C) override;

    inline void setNPar() override
    {
        nPar = 4;
    }

private:
    double B, E, G, Offset, m_Estart, m_Eend;
    bool isSubtracted;
};

#endif
