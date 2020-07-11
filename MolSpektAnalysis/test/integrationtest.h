//
// C++ Interface: IntegrationTest
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef INTEGRATIONTEST_H
#define INTEGRATIONTEST_H


class Potential;


class IntegrationTest
{
public:
    IntegrationTest(Potential* const potential);

    void runIntegrationTest(const double Rmin, const double Rmax, const double expectedFinalPrecision) const;

private:
    double testSingleIntegration(const double Rmin, const double Rmax, const int numPoints) const;

    Potential* const m_potential;
};

#endif // INTEGRATIONTEST_H
