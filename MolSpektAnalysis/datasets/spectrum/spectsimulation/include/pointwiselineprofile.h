//
// C++ Interface: PointwiseLineProfile
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POINTWISELINEPROFILE_H
#define POINTWISELINEPROFILE_H


#include "lineprofilepoint.h"


class QString;


class PointwiseLineProfile
{
public:
    PointwiseLineProfile();
    virtual ~PointwiseLineProfile();

    bool readData(const QString& Filename);

    inline int GetNumPoints() const
    {
        return m_NumPoints;
    }

    inline double GetEMax() const
    {
        return (m_NumPoints > 0 ? m_points[m_NumPoints - 1].Energy : 0.0);
    }

    inline double GetEnergy(const int i_N) const
    {
        return m_points[i_N].Energy;
    }

    inline double GetIntensity(const int i_N) const
    {
        return m_points[i_N].Intensity;
    }

private:
    LineProfilePoint* m_points;
    int m_NumPoints;
};

#endif
