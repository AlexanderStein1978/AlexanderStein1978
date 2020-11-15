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
#include "particle.h"


class CalculationTest : public ::testing::Test
{
public:
    CalculationTest()
        : Calc()
    {

    }

    ~CalculationTest()
    {

    }

    void SetUp()
    {

    }

    void TearDown()
    {

    }

protected:

    bool bindsParticleTo(Particle* P1, Particle* P2)
    {
        for (int n=0; n<4; ++n) if (P1->bound[n] == P2) return true;
        return false;
    }

    bool areParticlesBound(Particle* P1, Particle* P2)
    {
        return (bindsParticleTo(P1, P2) && bindsParticleTo(P2, P1));
    }

    bool areParticleBindingsCorrectlyInitialized(int index)
    {
        return (areParticlesBound(Calc.P + index, Calc.P + index - 2)
             && areParticlesBound(Calc.P + index, Calc.P + index + 2)
             && areParticlesBound(Calc.P + index, Calc.P + index - 2 * Calc.PXS)
             && areParticlesBound(Calc.P + index, Calc.P + index + 2 * Calc.PXS));
    }

    int getPXS()
    {
        return Calc.PXS;
    }

    int getPZS()
    {
        return Calc.PZS;
    }

    Calculation Calc;
};

TEST_F(CalculationTest, CheckParticleBindingInitialisation)
{
    for (int z=1; z < getPZS() - 1; ++z) for (int x=1; x < getPXS() - 1; ++x)
    {
        int index = 2 * (getPXS() * z + x);
        EXPECT_TRUE(areParticleBindingsCorrectlyInitialized(index)) << "index = " << index;
        EXPECT_TRUE(areParticleBindingsCorrectlyInitialized(index + 1)) << "index + 1 = " << (index + 1);
    }
}
