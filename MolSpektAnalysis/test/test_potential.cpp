//
// C++ Implementation: PotentialTest
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#include "integrationtest.h"
#include "gtest/gtest.h"
#include "potential.h"
#include "SplinePoint.h"


class PotentialTest : public ::testing::Test
{
public:
    PotentialTest()
        : pot(new Potential)
        , intTest(new IntegrationTest(pot))
    {

    }

    ~PotentialTest()
    {
        delete pot;
        delete intTest;
    }

    void SetUp()
    {
        ASSERT_TRUE(pot->readData("TestSpline.pot"));

        ASSERT_EQ(SplinePotential, pot->getPotType());
    }

    void TearDown()
    {

    }
    
    double getSplineSlope(const SplinePoint* const points, const int p, const double A, const double B)
    {
        double deltaX, dDeltaX = -1.0, Q, deltaXdsix, yssF1, yssF2;
        const double one = 1.0, three = 3.0, dsix = one / 6.0;
        dDeltaX = 1.0 / (deltaX = points[p].x - points[p-1].x);
        Q = (points[p].y - points[p-1].y) * dDeltaX;
        deltaXdsix = deltaX * dsix;
        yssF1 = deltaXdsix * points[p-1].yss;
        yssF2 = deltaXdsix * points[p].yss;
        return Q - (three * A * A - one) * yssF1 + (three * B * B - one) * yssF2;
    }

    double getLRPoint(const double R, const double UInf, const double *const LRC, const int *const PLRC, const int NLRC)
    {
        int p = PLRC[NLRC-1];
        double Res = LRC[NLRC - 1];
        for (int n = NLRC - 2; n>=0 ; n--)
        {
            while (p > PLRC[n])
            {
                p--;
                Res /= R;
            }
            Res += LRC[n];
        }
        Res *= pow(R, -1 * PLRC[0]);
        return UInf - Res;
    }

    double getLRSlope(const double R, const double *const LRC, const int *const PLRC, const int NLRC)
    {
        double Res, dR;
        const double one = 1.0;
        int p, i;
        for (i = NLRC - 2, p = PLRC[NLRC-1], Res = static_cast<double>(p) * LRC[NLRC - 1], dR = one / R; i>=0 ; --i)
        {
            while (p > PLRC[i])
            {
                --p;
                Res *= dR;
            }
            Res += static_cast<double>(PLRC[i]) * LRC[i];
        }
        return Res * pow(dR, one + static_cast<double>(PLRC[0]));
     }

protected:
    Potential* pot;
    IntegrationTest* intTest;
};

TEST_F(PotentialTest, CheckInnerWall)
{
    double RM(pot->getInnerConnectionRadius()), Rm(0.5 * RM);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 0.1);
}

TEST_F(PotentialTest, CheckSplineRegion)
{
    double Rm(pot->getInnerConnectionRadius()), RM(pot->getOuterConnectionRadius());
    ASSERT_LT(0.0, Rm);
    ASSERT_LT(Rm, RM);

    intTest->runIntegrationTest(Rm, RM, 1e-3);
}

TEST_F(PotentialTest, CheckLongRangePart)
{
    double Rm(pot->getOuterConnectionRadius()), RM(2.0 * Rm);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 1e-5);
}

TEST_F(PotentialTest, VerifyThatPotIsContinuouslyDifferentiable)
{
    double Rm(pot->getInnerConnectionRadius()), RM(pot->getOuterConnectionRadius()), *LRC, iA, iO, Exp;
    int NSplinePoints, NLRC, *pLRC;
    SplinePoint *points;
    pot->getSplinePotForWriting(NSplinePoints, points, NLRC, pLRC, LRC, iA, iO, Exp);
    pot->cdConnectSR();
    pot->cdConnectLR(pLRC[NLRC - 2]);
    
    EXPECT_DOUBLE_EQ(iO + iA * pow(Rm, -Exp), points[0].y);
    EXPECT_DOUBLE_EQ(-Exp * iA * pow(Rm, -Exp - 1.0), getSplineSlope(points, 1, 1.0, 0.0));
    
    EXPECT_DOUBLE_EQ(points[NSplinePoints - 1].y, getLRPoint(points[NSplinePoints - 1].x, pot->getUinf(), LRC, pLRC, NLRC));
    EXPECT_DOUBLE_EQ(getSplineSlope(points, NSplinePoints - 1, 0.0, 1.0), getLRSlope(points[NSplinePoints - 1].x, LRC, pLRC, NLRC));
}
