//
// C++ Interface: Lorentzian
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LORENTZIAN_H
#define LORENTZIAN_H


#include "lineprofile.h"


class Lorentzian : public LineProfile
{
public:
    Lorentzian();
	Lorentzian(double *x, double *y, double *Sig, int N);
    Lorentzian(const QString& data);

    void setData(double* x, double* y, double* Sig, int N) override;
    double GetPoint(double i_E) const override;
    double GetProfilePoint(double i_relE) const override;
    void Serialize(QTextStream& stream, const bool finish = true) const override;
    void getLineY(double *Ycalc) const override;

    inline double GetEcenter() const override
    {
        return E;
    }

    inline bool isValid() const override
    {
        return 0.0 != A && 0.0 != Gamma;
    }

    inline bool isWithSaturation() const override
    {
        return false;
    }
	
protected:
    void getDerivatives(double **deriv) override;
    bool getCalcY(double *YCalc) const override;
    void getPar(double *Par) const override;
    void setPar(double *Par) override;
    void updatePar(double *C) override;
	
private:
	double A, E, Gamma, AWF, EWF, GWF;
};

#endif
