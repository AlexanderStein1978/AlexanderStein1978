//
// C++ Interface: Gaussian line profile
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef GAUSSIAN_H
#define GAUSSIAN_H


#include "lineprofile.h"

class QString;


class Gaussian : public LineProfile
{
public:
    Gaussian();
    Gaussian(double *x, double *y, double *Sig, int N);
    Gaussian(const QString& data);
    virtual ~Gaussian();
    void GetValues(double& o_Intensity, double& o_CenterFreq, double& o_Width, double& o_Offset) const;
    void SetValues(const double Intensity, const double CenterFreq, const double Width, const double Offset);
    virtual double GetPoint(double i_E) const;
    virtual double GetProfilePoint(double i_relE) const;
    virtual void Serialize(QTextStream& stream, const bool finish = true) const;
    virtual void getLineY(double *Ycalc) const;

    inline double GetEcenter() const
    {
        return E;
    }

    inline bool isValid() const
    {
        return B != 0.0;
    }

    virtual inline bool isWithSaturation() const
    {
        return false;
    }

protected:
    bool getCalcYAndDerivatives(double *Ycalc, double **deriv) override;
    bool getCalcY(double *Ycalc) const override;
    void getPar(double *Par) const override;
    void setPar(double *Par) override;
    void updatePar(double *C) override;
    void initialize(const QStringList& data);

    inline void setNPar() override
    {
        nPar = 4;
    }

    double B, E, G;
};

#endif
