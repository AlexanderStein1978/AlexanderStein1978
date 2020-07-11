//
// C++ Implementation: PotentialTest
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

protected:
    Potential* pot;
    IntegrationTest* intTest;
};

TEST_F(PotentialTest, CheckInnerWall)
{
    double RM(pot->getInnerConnectionRadius()), Rm(0.5 * RM);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 1e-5);
}

TEST_F(PotentialTest, CheckSplineRegion)
{
    double Rm(pot->getInnerConnectionRadius()), RM(pot->getOuterConnectionRadius());
    ASSERT_LT(0.0, Rm);
    ASSERT_LT(Rm, RM);

    intTest->runIntegrationTest(Rm, RM, 1e-5);
}

TEST_F(PotentialTest, CheckLongRangePart)
{
    double Rm(pot->getOuterConnectionRadius()), RM(2.0 * Rm);
    ASSERT_LT(0.0, Rm);

    intTest->runIntegrationTest(Rm, RM, 1e-5);
}
