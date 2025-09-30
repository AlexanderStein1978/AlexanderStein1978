//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef GAUSSIANWITHSATURATION_H
#define GAUSSIANWITHSATURATION_H


#include "gaussian.h"

class Lorentzian;


class GaussianWithSaturation : public Gaussian
{
public:
    GaussianWithSaturation(double *x, double *y, double *Sig, int N);
    GaussianWithSaturation(const QString& data);
    GaussianWithSaturation(const Gaussian& other);
    GaussianWithSaturation(const Lorentzian& other);
    GaussianWithSaturation(const LineProfile& other);
    virtual ~GaussianWithSaturation();
    double GetPoint(double i_E) const override;
    void getLineY(double *Ycalc) const override;

    inline bool isWithSaturation() const override
    {
        return true;
    }

protected:
    bool getCalcYAndDerivatives(double *Ycalc, double **deriv) override;
    bool getCalcY(double *Ycalc) const override;
};

#endif // GAUSSIANWITHSATURATION_H
