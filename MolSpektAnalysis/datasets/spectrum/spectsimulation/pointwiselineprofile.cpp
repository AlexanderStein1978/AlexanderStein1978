//
// C++ Implementation: PointwiseLineProfile
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "pointwiselineprofile.h"

#include <QFile>
#include <QTextStream>


PointwiseLineProfile::PointwiseLineProfile()
    : m_points(NULL), m_NumPoints(0)
{
}

PointwiseLineProfile::~PointwiseLineProfile()
{
    if (m_points != NULL) delete[] m_points;
}

bool PointwiseLineProfile::readData(const QString &Filename)
{
    if (m_NumPoints > 0)
    {
        delete[] m_points;
        m_points = NULL;
        m_NumPoints = 0;
    }
    QFile file(Filename);
    bool Result = file.open(QIODevice::ReadOnly);
    if (Result)
    {
        QTextStream S(&file);
        QStringList rows;
        S.readLine();
        while (!S.atEnd()) rows.push_back(S.readLine());
        if (rows.size() > 0) m_points = new LineProfilePoint[rows.size()];
        for (QStringList::const_iterator it = rows.begin(); it != rows.end(); ++it)
        {
            QStringList row = it->split('\t');
            if (row.size() == 2)
            {
                m_points[m_NumPoints].Energy = row[0].toDouble();
                m_points[m_NumPoints++].Intensity = row[1].toDouble();
            }
        }
        if (m_NumPoints == 0)
        {
            if (m_points != NULL)
            {
                delete[] m_points;
                m_points = NULL;
            }
            Result = false;
        }
    }
    return Result;
}
