//
// C++ Implementation: LocalPerturbation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "LocalPerturbation.h"
#include "fit.h"

#include <QTextStream>


LocalPerturbation::LocalPerturbation()
    : m_bandConstants(), m_upOBCaBandConst(), m_H12(0.0, 0.0), m_Omega(0.0), m_IsoF(0.0), m_chisq(0.0), m_BeCurState(0.0),
      m_Center(0.0), m_data(), m_JStep(0), m_lastLevelUp(-1), m_lastLevelDown(-1), m_minJ(-1), m_maxJ(-1), m_firstUp(true), m_BeFixed(false)
{
}

LocalPerturbation::LocalPerturbation(const double i_IsoF, const double i_Omega, const std::vector<Js>& i_data, const int i_JStep,
                                     const double i_BeCurState)
    : m_bandConstants(), m_upOBCaBandConst(), m_H12(0.0, 0.0), m_Omega(i_Omega), m_IsoF(i_IsoF), m_chisq(0.0), m_BeCurState(i_BeCurState),
      m_Center(0.0), m_data(i_data), m_JStep(i_JStep), m_lastLevelUp(-1), m_lastLevelDown(-1), m_minJ(-1), m_maxJ(-1), m_firstUp(true), m_BeFixed(false)
{
    SetData(m_data);
    m_maxBadBest = 5;
}

LocalPerturbation::LocalPerturbation(const LocalPerturbation &i_toCopy)
    : QObject(), FitObject(i_toCopy), m_bandConstants(i_toCopy.m_bandConstants), m_upOBCaBandConst(i_toCopy.m_upOBCaBandConst), m_H12(i_toCopy.m_H12),
      m_Omega(i_toCopy.m_Omega), m_IsoF(i_toCopy.m_IsoF), m_chisq(i_toCopy.m_chisq), m_BeCurState(i_toCopy.m_BeCurState),
      m_Center(i_toCopy.m_Center), m_data(i_toCopy.m_data), m_JStep(i_toCopy.m_JStep), m_lastLevelUp(i_toCopy.m_lastLevelUp), m_lastLevelDown(i_toCopy.m_lastLevelDown),
      m_minJ(i_toCopy.m_minJ), m_maxJ(i_toCopy.m_maxJ), m_firstUp(i_toCopy.m_firstUp), m_BeFixed(i_toCopy.m_BeFixed)
{
    SetData(m_data);
    m_maxBadBest = 5;
}

LocalPerturbation::~LocalPerturbation()
{
}

LocalPerturbation& LocalPerturbation:: operator=(const LocalPerturbation& i_right)
{
    if (this != &i_right)
    {
        FitObject::operator=(i_right);
        m_bandConstants = i_right.m_bandConstants;
        m_upOBCaBandConst = i_right.m_upOBCaBandConst;
        m_H12 = i_right.m_H12;
        m_Omega = i_right.m_Omega;
        m_IsoF = i_right.m_IsoF;
        m_chisq = i_right.m_chisq;
        m_BeCurState = i_right.m_BeCurState;
        m_Center = i_right.m_Center;
        SetData(i_right.m_data);
        m_JStep = i_right.m_JStep;
        m_firstUp = i_right.m_firstUp;
        m_BeFixed = i_right.m_BeFixed;
    }
    return *this;
}

void LocalPerturbation::SetElStateSpecificValues(const double i_IsoF, const double i_Omega, const int i_JStep, const double i_BeCurState)
{
    m_IsoF = i_IsoF;
    m_Omega = i_Omega;
    m_JStep = i_JStep;
    m_BeCurState = i_BeCurState;
}

void LocalPerturbation::RemoveLastBandC()
{
    if (m_bandConstants.size() > 0) m_bandConstants.resize(m_bandConstants.size() - 1);
}

void LocalPerturbation::SetBandC(const int i_n, const double i_value, const double i_stdDev)
{
    if (i_n < static_cast<int>(m_bandConstants.size()))
    {
        m_bandConstants[i_n].value = i_value;
        m_bandConstants[i_n].stdDev = i_stdDev;
    }
    else
    {
        BandConst newC(i_value, i_stdDev);
        m_bandConstants.push_back(newC);
    }
}

void LocalPerturbation::CalcUpOBCaBandConst()
{
    bool fitOK = true;
    int NData = m_data.size();
    double JArray[NData], EArray[NData], uncArray[NData];
    for (int n=0; n < NData; ++n)
    {
        JArray[n] = m_data[n].J;
        EArray[n] = m_data[n].E_calc;
        uncArray[n] = m_data[n].unc;
    }
    for (int NBandC = 2; fitOK; ++NBandC)
    {
        double BandC[NBandC], err[NBandC];
        FitBandConstants(BandC, err, NBandC, m_IsoF, m_Omega, JArray, EArray, uncArray, NData);
        for (int n=0; n < NBandC; ++n) if (BandC[n] / err[n] < 10.0) fitOK = false;
        if (fitOK || NBandC == 2)
        {
            m_upOBCaBandConst.clear();
            for (int n=0; n < NBandC; ++n)
            {
                BandConst newC(BandC[n], err[n]);
                m_upOBCaBandConst.push_back(newC);
            }
        }
    }
}

double LocalPerturbation::GetPointFromBandAndJFactor(const std::vector<BandConst>& i_bandConst, const double i_JFactor)
{
    double Res = i_bandConst.back().value;
    for (std::vector<BandConst>::const_reverse_iterator it = i_bandConst.rbegin() + 1; it != i_bandConst.rend(); ++it)
    {
        Res *= i_JFactor;
        Res += it->value;
    }
    return Res;
}

double LocalPerturbation::GetPointPert(const bool i_up, const double i_JValue) const
{
    double F = m_IsoF * (i_JValue * (i_JValue + 1.0) - m_Omega * m_Omega), E1 = GetPointFromBandAndJFactor(m_upOBCaBandConst, F), E2 = GetPointFromBandAndJFactor(m_bandConstants, F), halfDiff = 0.5 * (E1 - E2);
    double root = sqrt(halfDiff * halfDiff + m_H12.value * m_H12.value);
    return 0.5 * (E1 + E2) + (i_up ? root : - root);
}

double LocalPerturbation::GetEnergyDiff(const double i_JValue) const
{
    double F = m_IsoF * (i_JValue * (i_JValue + 1.0) - m_Omega * m_Omega), E1 = GetPointFromBandAndJFactor(m_upOBCaBandConst, F), E2 = GetPointFromBandAndJFactor(m_bandConstants, F), halfDiff = 0.5 * (E1 - E2);
    double root = sqrt(halfDiff * halfDiff + m_H12.value * m_H12.value);
    return 0.5 * (E1 - E2) + (E2 > E1 ? root : - root);
}

double LocalPerturbation::GetFitSigma() const
{
    return sqrt(m_chisq / (m_data.size() - m_bandConstants.size() - 1));
}

bool LocalPerturbation::getCalcYAndDerivatives(double *o_ycalc, double **o_deriv)
{
    double Hsq = m_H12.value * m_H12.value, JtJp1 = m_IsoF * (m_data[0].J * (m_data[0].J + 1.0) - m_Omega * m_Omega);
    double E2_calc, DeltaE_calc, lRoot, throughlRoot, Ecalc_u, Ecalc_d;
    bool RetVal = true;
    for (unsigned int n=0; n < m_data.size(); ++n)
    {
        if (n==0 || m_data[n-1].J < m_data[n].J)
        {
            if (n>0) for (int j = static_cast<int>(m_data[n-1].J) + 1; j <= static_cast<int>(m_data[n].J); ++j) JtJp1 += m_IsoF * static_cast<double>(2*j);
            std::vector<BandConst>::reverse_iterator it = m_bandConstants.rbegin();
            E2_calc = it->value;
            for (++it; it != m_bandConstants.rend(); ++it)
            {
                E2_calc *= JtJp1;
                E2_calc += it->value;
            }
            DeltaE_calc = m_data[n].E_calc - E2_calc;
            lRoot = sqrt(0.25 * DeltaE_calc * DeltaE_calc + Hsq);
            throughlRoot = 1.0 / lRoot;
            double E_av = 0.5 * (m_data[n].E_calc + E2_calc);
            Ecalc_u = E_av + lRoot;
            Ecalc_d = E_av - lRoot;
        }
        bool upper = abs(Ecalc_u - m_data[n].E_obs) < abs(Ecalc_d - m_data[n].E_obs);
        if (upper != m_data[n].isUp) RetVal = false;
        o_ycalc[n] = (upper ? Ecalc_u : Ecalc_d);
        o_deriv[n][0] = m_H12.value * (upper ? throughlRoot : -throughlRoot);
        o_deriv[n][1] = 0.5 + (upper ? -0.25 : 0.25) * DeltaE_calc * throughlRoot;
        if (!m_BeFixed) for (unsigned int i=1; i < m_bandConstants.size(); ++i) o_deriv[n][i+1] = o_deriv[n][i] * JtJp1;
        else if (m_bandConstants.size() > 2)
        {
            o_deriv[n][2] = o_deriv[n][1] * JtJp1 * JtJp1;
            for (unsigned int i=3; i < m_bandConstants.size(); ++i) o_deriv[n][i] = o_deriv[n][i-1] * JtJp1;
        }
        if (0 != m_DebugLogStream)
        {
            *m_DebugLogStream << QString("H=%1\tT=%2\tB=%3\tNBandC=%4\tJ=%5\tEcalc=%6\t_Eobs=%7\n").arg(m_H12.value, 0, 'f', 4)
                                 .arg(m_bandConstants[0].value, 0, 'f', 4).arg(m_bandConstants[1].value, 0, 'f', 6)
                                 .arg(m_bandConstants.size()).arg(m_data[n].J).arg(o_ycalc[n], 0, 'f', 4).arg(m_data[n].E_obs, 0, 'f', 4);
        }
    }
    return (RetVal && m_firstUp == (m_bandConstants[1].value > m_BeCurState));
}

void LocalPerturbation::getPar(double *o_par) const
{
    o_par[0] = m_H12.value;
    if (m_BeFixed)
    {
        o_par[1] = m_bandConstants[0].value;
        for (unsigned int n=2; n < m_bandConstants.size(); ++n) o_par[n] = m_bandConstants[n].value;
    }
    else for (unsigned int n=0; n < m_bandConstants.size(); ++n) o_par[n+1] = m_bandConstants[n].value;
}

void LocalPerturbation::setPar(double *i_par)
{
    m_H12.value = i_par[0];
    if (m_BeFixed)
    {
        m_bandConstants[0].value = i_par[1];
        for (unsigned int n=2; n < m_bandConstants.size(); ++n) m_bandConstants[n].value = i_par[n];
    }
    else for (unsigned int n=0; n < m_bandConstants.size(); ++n) m_bandConstants[n].value = i_par[n+1];
}

void LocalPerturbation::updatePar(double *i_par)
{
    m_H12.value += i_par[0];
    if (m_BeFixed)
    {
        m_bandConstants[0].value += i_par[1];
        for (unsigned int n=2; n < m_bandConstants.size(); ++n) m_bandConstants[n].value += i_par[n];
    }
    else for (unsigned int n=0; n < m_bandConstants.size(); ++n) m_bandConstants[n].value += i_par[n+1];
}

void LocalPerturbation::setNPar()
{
    nPar = m_bandConstants.size() + (m_BeFixed ? 0 : 1);
}

void LocalPerturbation::SetData(const std::vector<Js>& i_data)
{
    m_data = i_data;
    int N = i_data.size();
    double *x = new double[N], *y = new double[N], *sig = new double[N];
    for (int n=0; n<N; ++n)
    {
        x[n] = i_data[n].J;
        y[n] = i_data[n].E_obs;
        sig[n] = 1.0 / (i_data[n].unc * i_data[n].unc);
    }
    setData(x, y, sig, N);
}

void LocalPerturbation::writeData(QTextStream *i_stream) const
{
    *i_stream << "Begin LocalPerturbation\n";
    *i_stream << "Min J: " << QString::number(m_minJ >= 0 ? m_minJ : (m_data.size() > 0 ? m_data[0].J : -1)) << '\n';
    *i_stream << "Max J: " << QString::number(m_maxJ >= 0 ? m_maxJ : (m_data.size() > 0 ? m_data[m_data.size() - 1].J : -1)) << '\n';
    *i_stream << "Center J: " << QString::number(m_Center, 'f', 1) << '\n';
    *i_stream << "H_12: " << QString::number(m_H12.value, 'g', 12) << '(' << QString::number(m_H12.stdDev, 'g', 6) << ')' << '\n';
    *i_stream << "Number of band constants: " << QString::number(m_bandConstants.size()) << '\n';
    for (unsigned int n=0; n < m_bandConstants.size(); ++n)
        *i_stream << QString::number(n) << ": " << QString::number(m_bandConstants[n].value, 'g', 12)
                  << '(' << QString::number(m_bandConstants[n].stdDev, 'g', 6) << ')' << '\n';
    *i_stream << "ChiSq: " << QString::number(m_chisq, 'f', 3) << '\n';
    *i_stream << "End LocalPerturbation\n";
}

bool LocalPerturbation::readData(QTextStream *i_stream)
{
    bool EndFound = false, Result = true;
    int m, n, o, i, N = 0, Nf = 0;
    while (!i_stream->atEnd())
    {
        QString Buffer = i_stream->readLine();
        if (Buffer.indexOf("Min J:", 0, Qt::CaseInsensitive) >= 0) m_minJ = Buffer.right(Buffer.length() - 6).toInt();
        else if (Buffer.indexOf("Max J:", 0, Qt::CaseInsensitive) >= 0) m_maxJ = Buffer.right(Buffer.length() - 6).toInt();
        else if (Buffer.indexOf("Center J:", 0, Qt::CaseInsensitive) >= 0) m_Center = Buffer.right(Buffer.length() - 9).toDouble();
        else if ((m = Buffer.indexOf("H_12:", 0, Qt::CaseInsensitive)) >= 0)
        {
            n = Buffer.indexOf('(', m);
            o = Buffer.indexOf(')', n);
            if (n==-1) n = Buffer.length();
            if (n - m > 5)
            {
                 if ((m_H12.value = Buffer.mid(m + 5, n - m - 5).toDouble()) == 0.0) Result = false;
            }
            else Result = false;
            m_H12.stdDev = (o - n > 1 ? Buffer.mid(n+1, o-n-1).toDouble() : 0.0);
        }
        else if (Buffer.indexOf("Number of band constants:", 0, Qt::CaseInsensitive) >= 0)
        {
            N = Buffer.right(Buffer.length() - 25).toInt();
            m_bandConstants.resize(N);
            for (std::vector<BandConst>::iterator it = m_bandConstants.begin(); it != m_bandConstants.end(); ++it)
            {
                it->value = 0.0;
                it->stdDev = 0.0;
            }
        }
        else if (Buffer.indexOf("ChiSq:", 0, Qt::CaseInsensitive) >= 0) m_chisq = Buffer.right(Buffer.length() - 6).toDouble();
        else if (Buffer.indexOf("End LocalPerturbation", 0, Qt::CaseInsensitive) >= 0)
        {
            EndFound = true;
            break;
        }
        else if ((m = Buffer.indexOf(':')) >= 1)
        {
            n = Buffer.indexOf('(', m);
            o = Buffer.indexOf(')', n);
            i = Buffer.left(m).toInt();
            if (n==-1) n = Buffer.length();
            if (i < 0 || i >= N || m_bandConstants[i].value != 0.0 || n - m <= 1) Result = false;
            else
            {
                m_bandConstants[i].value = Buffer.mid(m+1, n-m-1).toDouble();
                ++Nf;
                if (o - n > 1) m_bandConstants[i].stdDev = Buffer.mid(n+1, o-n-1).toDouble();
            }
        }
    }
    if (!EndFound || Nf < N) Result = false;
    return Result;
}

double LocalPerturbation::LevenbergMarquardt(int MaxIt, double MinImp)
{
    m_chisq = FitObject::LevenbergMarquardt(MaxIt, MinImp);
    emit Changed();
    return m_chisq;
}

void LocalPerturbation::SetLastLevels(const int i_lastUp, const int i_lastDown)
{
    m_lastLevelUp = i_lastUp;
    m_lastLevelDown = i_lastDown;
}
