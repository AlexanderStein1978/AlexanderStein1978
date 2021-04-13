//
// C++ Interface: Lorentzian with saturation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2021 - 2021
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LORENTZIANWITHSATURATION_H
#define LORENTZIANWITHSATURATION_H


#include "lorentzian.h"

class Gaussian;


class LorentzianWithSaturation : public Lorentzian
{
public:
    LorentzianWithSaturation();
    LorentzianWithSaturation(double *x, double *y, double *Sig, int N);
    LorentzianWithSaturation(const QString& data);
    LorentzianWithSaturation(const Gaussian& other);
    LorentzianWithSaturation(const Lorentzian& other);
    LorentzianWithSaturation(const LineProfile& other);

    double GetPoint(double i_E) const override;
    void getLineY(double *Ycalc) const override;

    inline bool isWithSaturation() const override
    {
        return true;
    }

protected:
    bool getCalcYAndDerivatives(double *Ycalc, double **deriv) override;
    void getDerivatives(double **deriv) override;
    bool getCalcY(double *Ycalc) const override;
};

#endif // LORENTZIANWITHSATURATION_H
