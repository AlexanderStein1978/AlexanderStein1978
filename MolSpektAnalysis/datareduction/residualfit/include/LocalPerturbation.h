//
// C++ Interface: LocalPerturbation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#ifndef LOCALPERTURBATION_H
#define LOCALPERTURBATION_H


#include "fitobject.h"
#include "BandConst.h"
#include "Js.h"

#include <QObject>


class LocalPerturbation : public QObject, public FitObject
{
    Q_OBJECT

public:
    LocalPerturbation();
    LocalPerturbation(const double i_IsoF, const double i_Omega, const std::vector<Js>& i_data, const int i_JStep, const double i_BeCurState);
    LocalPerturbation(const LocalPerturbation& i_toCopy);
    LocalPerturbation& operator=(const LocalPerturbation& i_right);
    virtual ~LocalPerturbation();
    void SetElStateSpecificValues(const double i_IsoF, const double i_Omega, const int i_JStep, const double i_BeCurState);
    void SetBandC(const int i_n, const double i_value, const double i_stdev);
    void RemoveLastBandC();
    double GetFitSigma() const;
    double GetEnergyDiff(const double i_JValue) const;
    double GetPointPert(bool i_up, const double i_JFactor) const;
    void SetData(const std::vector<Js>& i_data);
    double LevenbergMarquardt(int MaxIt, double MinImp);
    void SetLastLevels(const int i_lastUp, const int i_lastDown);
    void writeData(QTextStream* i_stream) const;
    bool readData(QTextStream* i_stream);

    inline double GetPointUpOBCa(const double i_JValue)
    {
        if (m_upOBCaBandConst.empty()) CalcUpOBCaBandConst();
        return GetPointFromBandAndJ(m_upOBCaBandConst, i_JValue);
    }

    inline double GetPointCalcBand(const double i_JValue) const
    {
        return GetPointFromBandAndJ(m_bandConstants, i_JValue);
    }

    double GetPointLowPert(const double i_JValue) const
    {
        return GetPointPert(false, i_JValue);
    }

    double GetPointUpPert(const double i_JValue) const
    {
        return GetPointPert(true, i_JValue);
    }

    inline void GetH12(double& o_H12, double& o_stdev) const
    {
        o_H12 = m_H12.value;
        o_stdev = m_H12.stdDev;
    }

    inline void SetH12(const double i_H12, const double i_stdev)
    {
        m_H12.value = i_H12;
        m_H12.stdDev = i_stdev;
    }

    inline void SetH12Uncertainty(const double i_uncertainty)
    {
        m_H12.stdDev = i_uncertainty;
    }

    inline int GetNBandC() const
    {
        return m_bandConstants.size();
    }

    inline void GetBandC(const int i_n, double& o_value, double& o_stdDev) const
    {
        o_value = m_bandConstants[i_n].value;
        o_stdDev = m_bandConstants[i_n].stdDev;
    }

    inline void SetBandCUncertainty(const int i_n, const double i_uncertainty)
    {
        m_bandConstants[i_n].stdDev = i_uncertainty;
    }

    inline int GetNData() const
    {
        return m_data.size();
    }

    inline double GetDataJValue(const int i_n) const
    {
        return static_cast<double>(m_data[i_n].J);
    }

    inline double GetDataEobs(const int i_n) const
    {
        return m_data[i_n].E_obs;
    }

    inline double GetDataECalc(const int i_n) const
    {
        return m_data[i_n].E_calc;
    }

    inline double GetDataUnc(const int i_n) const
    {
        return m_data[i_n].unc;
    }

    inline int GetJStep() const
    {
        return m_JStep;
    }

    inline double GetOmega() const
    {
        return m_Omega;
    }

    inline double GetIsoF() const
    {
        return m_IsoF;
    }

    inline void SetChisq(const double i_chisq)
    {
        m_chisq = i_chisq;
        m_upOBCaBandConst.clear();
    }

    inline double GetChisq() const
    {
        return m_chisq;
    }

    inline double GetBeCurState() const
    {
        return m_BeCurState;
    }

    inline double GetCenter() const
    {
        return m_Center;
    }

    inline void SetCenter(const double i_center)
    {
        m_Center = i_center;
    }

    inline void SetFirstUp(const bool i_firstUp)
    {
        m_firstUp = i_firstUp;
    }

    inline int GetLastLevelUp() const
    {
        return m_lastLevelUp;
    }

    inline int GetLastLevelDown() const
    {
        return m_lastLevelDown;
    }

    inline void SetDataIsUp(const int i_n, const bool i_isUp)
    {
        m_data[i_n].isUp = i_isUp;
    }

    inline bool IsDataUp(const int i_n) const
    {
        return m_data[i_n].isUp;
    }

    inline void SetBeFixed(const bool i_fixed)
    {
        m_BeFixed = i_fixed;
    }

    inline int GetMinJ() const
    {
        return m_minJ;
    }

    inline int GetMaxJ() const
    {
        return m_maxJ;
    }

    inline std::vector<Js>& GetFitData()
    {
        return m_data;
    }

protected:
    virtual bool getCalcYAndDerivatives(double *o_ycalc, double **deriv) override;
    virtual void getPar(double *Par) const override;
    virtual void setPar(double *Par) override;
    virtual void setNPar() override;
    virtual void updatePar(double *C) override;

signals:
    void Changed();

private:
    void CalcUpOBCaBandConst();
    static double GetPointFromBandAndJFactor(const std::vector<BandConst>& i_bandConst, const double i_JFactor);

    inline double GetPointFromBandAndJ(const std::vector<BandConst>& i_bandConst, const double i_JValue) const
    {
        return GetPointFromBandAndJFactor(i_bandConst, m_IsoF * (i_JValue * (i_JValue + 1.0) - m_Omega * m_Omega));
    }

    std::vector<BandConst> m_bandConstants, m_upOBCaBandConst;
    BandConst m_H12;
    double m_Omega, m_IsoF, m_chisq, m_BeCurState, m_Center;
    std::vector<Js> m_data;
    int m_JStep, m_lastLevelUp, m_lastLevelDown, m_minJ, m_maxJ;
    bool m_firstUp, m_BeFixed;
};

#endif
