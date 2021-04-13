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

class Gaussian;


class Lorentzian : public LineProfile
{
public:
    Lorentzian();
	Lorentzian(double *x, double *y, double *Sig, int N);
    Lorentzian(const QString& data);
    Lorentzian(const Gaussian& other);
    Lorentzian(const Lorentzian& other);
    Lorentzian(const LineProfile& other);

    void setData(double* x, double* y, double* Sig, int N) override;
    double GetPoint(double i_E) const override;
    double GetProfilePoint(double i_relE) const override;
    void Serialize(QTextStream& stream, const bool finish = true) const override;
    void getLineY(double *Ycalc) const override;
    double GetWidth() const;
    void GetValues(double& o_Intensity, double& o_CenterFreq, double& o_Width, double& o_Offset) const override;
    void SetValues(const double Intensity, const double CenterFreq, const double Width, const double Offset) override;

    inline double GetECenter() const override
    {
        return E;
    }

    inline double GetAmplitude() const
    {
        return A;
    }

    inline bool isValid() const override
    {
        return 0.0 != A && 0.0 != Gamma;
    }

    inline bool isWithSaturation() const override
    {
        return false;
    }

    inline LineProfileType getType() const override
    {
        return LorentzianType;
    }

protected:
    void getDerivatives(double **deriv) override;
    bool getCalcY(double *YCalc) const override;
    void getPar(double *Par) const override;
    void setPar(double *Par) override;
    void updatePar(double *C) override;

    inline void setNPar() override
    {
        nPar = 4;
    }

private:
    void InitializeFromLineProfile(const LineProfile& other);

	double A, E, Gamma, AWF, EWF, GWF;    
};

#endif
