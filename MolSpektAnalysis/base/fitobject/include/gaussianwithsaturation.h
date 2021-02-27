//
// C++ Interface: Gaussian line profile with saturation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2014 - 2021
//
// Copyright: See README file that comes with this source code
//
//

#ifndef GAUSSIANWITHSATURATION_H
#define GAUSSIANWITHSATURATION_H


#include "gaussian.h"


class GaussianWithSaturation : public Gaussian
{
public:
    GaussianWithSaturation(double *x, double *y, double *Sig, int N);
    GaussianWithSaturation(const QString& data);
    virtual ~GaussianWithSaturation();
    void GetValues(double& o_Intensity, double& o_SatIntensity, double& o_SatBorder, double& o_CenterFreq, double& o_Width, double& o_Offset) const;
    void SetValues(const double Intensity, const double SatIntensity, const double SatBorder, const double CenterFreq, const double Width, const double Offset);
    double GetPoint(double i_E) const override;
    double GetProfilePoint(double i_relE) const override;
    void Serialize(QTextStream& stream, const bool finish = true) const override;
    void getLineY(double *Ycalc) const override;
    void setData(double *x, double *y, double *Sig, int N) override;

    inline bool isWithSaturation() const override
    {
        return true;
    }

protected:
    bool getCalcYAndDerivatives(double *Ycalc, double **deriv) override;
    bool getCalcY(double *Ycalc) const override;
    void getPar(double *Par) override;
    void setPar(double *Par) override;
    void updatePar(double *C) override;

    inline void setNPar() override
    {
        nPar = 5;
    }

private:
    enum WhichLine{LineBestfit, LineGaussian, LineSaturated};
    
    int getNBorderByI(const double bY, const int nStart, const int sDirection) const;
    int findBorder(const int start, const int direction) const;
    double GetPoint(double i_E, const WhichLine l) const;
    void updateBorders();

    double SB, IBo;
    int iBn, oBn;
};

#endif // GAUSSIANWITHSATURATION_H
