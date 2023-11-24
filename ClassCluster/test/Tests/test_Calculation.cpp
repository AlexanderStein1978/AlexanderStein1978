//
// C++ Implementation: CalculationTest
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "gtest/gtest.h"
#include "Calculation.h"
#include "CalculationTestHelper.h"
#include "particle.h"
#include "potstruct.h"
#include "potential.h"

#include <algorithm>


class CalculationTest : public ::testing::Test
{
public:
    CalculationTest()
        : Calc(nullptr)
        , cth(nullptr)
    {

    }

    ~CalculationTest()
    {
        if (nullptr != Calc) delete Calc;
    }

    void SetUp()
    {
        PotStruct struc[Calculation::NumPot];
        QString dataDir(DATA_DIRECTORY);
        struc[Calculation::ClosestTwo].InitAsOwner(dataDir + "/ClosestTwo.pot");
        struc[Calculation::NextTwo].InitAsOwner(dataDir + "/NextTwo.pot");
        struc[Calculation::Remaining].InitAsOwner(dataDir + "/Remaining.pot");
        struc[Calculation::Angular].InitAsOwner(dataDir + "/Angular.pot");
        Calc = new Calculation(struc);
        cth = new CalculationTestHelper(Calc);
    }

    void TearDown()
    {

    }

protected:

    Calculation* Calc;
    CalculationTestHelper* cth;
};


TEST_F(CalculationTest, CheckParticleBindingInitialisation)
{
    bool even = true;
    int start = 0;
    for (int z=1; z < cth->getPZS() - 1; ++z)
    {
        even = (even ? false : true);
        start += (even ? 38 : 40);
        for (int x=2; x < (even ? cth->getPXS() - 2 : cth->getPXS() - 4); ++x)
        {
            int index = start + x;
            EXPECT_TRUE(cth->areParticleBindingsCorrectlyInitialized(index)) << "index = " << index;
        }
    }
    for (int n=0; n < cth->getNumParticles(); ++n) EXPECT_EQ(cth->getParticleNB(n), cth->getParticleMNB(n)) << "index = " << n;
}

TEST_F(CalculationTest, CheckParticleBindingUpdates)
{
    int indexRowLength = 2 * cth->getPXS();
    cth->swapParticlePositions(indexRowLength + 2, indexRowLength + 4);
    cth->swapParticlePositions(10 * indexRowLength + 2, 11 * indexRowLength + 2);
    cth->updateBindings();

    for (int i=0; i < cth->getNumParticles(); ++i)
    {
        EXPECT_TRUE(cth->isNoBindingDoubled(i));
        
        for (int j=0; j < cth->getParticleNB(i); ++j)
        {
            int k = cth->getBoundParticleIndex(i, j);
            if (k<0) break;
            
            EXPECT_TRUE(cth->areParticlesBound(i, k));
            EXPECT_LT(cth->getDist(i, k), 6.0);
        }
    }
}

TEST_F(CalculationTest, CalcEndpointsOfEnergyDefinitionAxis)
{
    Vector end1, end2, direction(1.0, 1.0, 1.0);
    Calc->CalcEndpointsOfEnergyDefinitionAxis(210, direction, end1, end2);
    EXPECT_NEAR(12.9474, end1.X(), 1e-4);
    EXPECT_NEAR(16.9474, end1.Y(), 1e-4);
    EXPECT_EQ(0.0, end1.Z());
    EXPECT_NEAR(76.0, end2.X(), 1e-4);
    EXPECT_EQ(80.0, end2.Y());
    EXPECT_NEAR(63.0526, end2.Z(), 1e-4);
}

TEST_F(CalculationTest, SetEnergy_equalDistrib)
{
    EXPECT_EQ(0.0, Calc->getKineticEnergy());
    Calc->setEnergy(1e6);
    EXPECT_NEAR(1e6, Calc->getKineticEnergy(), 1e-7);
    Calc->setEnergy(1e6);
    EXPECT_NEAR(2e6, Calc->getKineticEnergy(), 1e-7);
    Calc->setEnergy(-1e6);
    EXPECT_NEAR(1e6, Calc->getKineticEnergy(), 1e-7);
    Calc->setEnergy(-1e6 - 1);
    EXPECT_NEAR(0.0, Calc->getKineticEnergy(), 1e-7);
}

TEST_F(CalculationTest, SetEnergy_randomDistrib)
{
    double kinEnergy = cth->createRandomSpeedDistribution();
    Calc->setEnergy(-0.7 * kinEnergy);
    EXPECT_NEAR(0.3 * kinEnergy, Calc->getKineticEnergy(), 1e-5);
}

TEST_F(CalculationTest, Triple)
{
    Vector R[3] = {Vector(32.0, 40.0, 40.0), Vector(40.0, 41.0, 40.0), Vector(48.0, 40.0, 40.0)}, v[3];
    int MNB[3] = {1, 2, 1};
    cth->setParticles(3, R, v, MNB);
    cth->addParticleBinding(0, 1);
    cth->addParticleBinding(1, 2);
    double EB = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CB = cth->getCenterOfMass(), vB = cth->getAverageV(), LB = cth->getAngularMomentum(CB);
    cth->run(3);
    double EA = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CA = cth->getCenterOfMass(), vA = cth->getAverageV(), LA = cth->getAngularMomentum(CA);
    EXPECT_NEAR(EA, EB, 1e-2);
    EXPECT_NEAR(CA.X(), CB.X(), 1e-2);
    EXPECT_NEAR(CA.Y(), CB.Y(), 1e-2);
    EXPECT_NEAR(CA.Z(), CB.Z(), 1e-2);
    EXPECT_NEAR(vA.X(), vB.X(), 1e-2);
    EXPECT_NEAR(vA.Y(), vB.Y(), 1e-2);
    EXPECT_NEAR(vA.Z(), vB.Z(), 1e-2);
    EXPECT_NEAR(LA.X(), LB.X(), 1e-2);
    EXPECT_NEAR(LA.Y(), LB.Y(), 1e-2);
    EXPECT_NEAR(LA.Z(), LB.Z(), 1e-2);
}

TEST_F(CalculationTest, Triple_FT)
{
    Vector R[3] = {Vector(32.000967, 40.124414, 40.0), Vector(40.0, 40.0, 40.0), Vector(47.999033, 40.124414, 40.0)}, v[3];
    int MNB[3] = {1, 2, 1};
    cth->setParticles(3, R, v, MNB);
    cth->addParticleBinding(0, 1);
    cth->addParticleBinding(1, 2);
    double EB = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CB = cth->getCenterOfMass(), vB = cth->getAverageV(), LB = cth->getAngularMomentum(CB);
    cth->run(3);
    double EA = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CA = cth->getCenterOfMass(), vA = cth->getAverageV(), LA = cth->getAngularMomentum(CA);
    EXPECT_NEAR(EA, EB, 1e-2);
    EXPECT_NEAR(CA.X(), CB.X(), 1e-2);
    EXPECT_NEAR(CA.Y(), CB.Y(), 1e-2);
    EXPECT_NEAR(CA.Z(), CB.Z(), 1e-2);
    EXPECT_NEAR(vA.X(), vB.X(), 1e-2);
    EXPECT_NEAR(vA.Y(), vB.Y(), 1e-2);
    EXPECT_NEAR(vA.Z(), vB.Z(), 1e-2);
    EXPECT_NEAR(LA.X(), LB.X(), 1e-2);
    EXPECT_NEAR(LA.Y(), LB.Y(), 1e-2);
    EXPECT_NEAR(LA.Z(), LB.Z(), 1e-2);
}
