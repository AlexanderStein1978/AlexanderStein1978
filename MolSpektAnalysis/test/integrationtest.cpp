//
// C++ Implementation: IntegrationTest
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "integrationtest.h"
#include "gtest/gtest.h"
#include "potential.h"

#include <cmath>


IntegrationTest::IntegrationTest(Potential* const potential)
    : m_potential(potential)
{

}

double IntegrationTest::testSingleIntegration(const double Rmin, const double Rmax, const int numPoints) const
{
    const double* const points = m_potential->get_dVdR(Rmin, Rmax, numPoints);
    double sum = 0.0;
    for (int n=0; n < numPoints; ++n) sum += points[n];
    EXPECT_FALSE(isnan(sum));
    delete[] points;
    return m_potential->getPoint(Rmax) - m_potential->getPoint(Rmin) - sum * (Rmax - Rmin) / numPoints;
}

void IntegrationTest::runIntegrationTest(const double Rmin, const double Rmax, const double expectedFinalPrecision) const
{
    int numPoints(20);
    for (double deviation = testSingleIntegration(Rmin, Rmax, numPoints / 2); abs(deviation) > expectedFinalPrecision; numPoints *= 2)
    {
        double lastDeviation = deviation;
        deviation = testSingleIntegration(Rmin, Rmax, numPoints);
        printf("numPoints = %d, deviation = %f\n", numPoints, deviation);
        ASSERT_GT(0.5 * abs(lastDeviation), abs(deviation));
    }
}

