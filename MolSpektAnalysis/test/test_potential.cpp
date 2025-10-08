//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
    }

    void TearDown()
    {
    }

    void LoadPotential(const QString fileName)
    {
        ASSERT_TRUE(pot->readData(fileName));

        ASSERT_EQ(SplinePotential, pot->getPotType());
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

    void verifyContinuousDifferentabilityOfSplinePotential()
    {
        double *LRC, iA, iO, Exp, prec(1e-12);
        int NSplinePoints, NLRC, *pLRC;
        SplinePoint *points;
        pot->cdConnectSR();
        pot->getSplinePotForWriting(NSplinePoints, points, NLRC, pLRC, LRC, iA, iO, Exp);
        pot->cdConnectLR(pLRC[NLRC - 2]);

        double ISL = pot->getSplineSlope(1, 1.0, 0.0), OSL = pot->getSplineSlope(NSplinePoints - 1, 0.0, 1.0);
        EXPECT_NEAR(iO + iA * pow(points[0].x, -Exp), points[0].y, abs(prec * points[0].y));
        EXPECT_NEAR(-Exp * iA * pow(points[0].x, -Exp - 1.0), ISL, abs(prec * ISL));

        EXPECT_NEAR(points[NSplinePoints - 1].y, getLRPoint(points[NSplinePoints - 1].x, pot->getUinf(), LRC, pLRC, NLRC), abs(prec * points[NSplinePoints - 1].y));
        EXPECT_NEAR(OSL, getLRSlope(points[NSplinePoints - 1].x, LRC, pLRC, NLRC), abs(prec * OSL));
    }

protected:
    Potential* pot;
    IntegrationTest* intTest;
};

TEST_F(PotentialTest, CheckInnerWall)
{
    LoadPotential("TestSpline.pot");
    double RM(pot->getInnerConnectionRadius()), Rm(0.5 * RM);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 0.1);
}

TEST_F(PotentialTest, CheckSplineRegion)
{
    LoadPotential("TestSpline.pot");
    double Rm(pot->getInnerConnectionRadius()), RM(pot->getOuterConnectionRadius());
    ASSERT_LT(0.0, Rm);
    ASSERT_LT(Rm, RM);

    intTest->runIntegrationTest(Rm, RM, 1e-3);
}

TEST_F(PotentialTest, CheckLongRangePart)
{
    LoadPotential("TestSpline.pot");
    double Rm(pot->getOuterConnectionRadius()), RM(2.0 * Rm);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 1e-5);
}

TEST_F(PotentialTest, VerifyThatPotIsContinuouslyDifferentiable)
{
    LoadPotential("TestSpline.pot");
    verifyContinuousDifferentabilityOfSplinePotential();
}

TEST_F(PotentialTest, VerifyThatPot2IsContinuouslyDifferentiable)
{
    LoadPotential("TestSpline2.pot");
    verifyContinuousDifferentabilityOfSplinePotential();
}
