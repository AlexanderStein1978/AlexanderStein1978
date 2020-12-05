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

#include <algorithm>


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

    bool bindsParticleTo(Particle* P1, Particle* P2) const
    {
        for (int n=0; n<4; ++n) if (P1->bound[n] == P2) return true;
        return false;
    }

    bool areParticlesBound(Particle* P1, Particle* P2) const
    {
        return (bindsParticleTo(P1, P2) && bindsParticleTo(P2, P1));
    }

    bool areParticlesBound(const int index1, const int index2) const
    {
        return areParticlesBound(Calc.P + index1, Calc.P + index2);
    }

    bool areParticleBindingsCorrectlyInitialized(int index)
    {
        return (areParticlesBound(Calc.P + index, Calc.P + index - 2)
             && areParticlesBound(Calc.P + index, Calc.P + index + 2)
             && areParticlesBound(Calc.P + index, Calc.P + index - 2 * Calc.PXS)
             && areParticlesBound(Calc.P + index, Calc.P + index + 2 * Calc.PXS));
    }

    void swapParticlePositions(const int index1, const int index2)
    {
        std::swap(Calc.P[index1].X, Calc.P[index2].X);
        std::swap(Calc.P[index1].Y, Calc.P[index2].Y);
        std::swap(Calc.P[index1].Z, Calc.P[index2].Z);
    }

    void updateBindings()
    {
        Calc.updateBindings();
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

TEST_F(CalculationTest, CheckParticleBindingUpdates)
{
    int indexRowLength = 2 * getPXS();
    swapParticlePositions(indexRowLength + 2, indexRowLength + 4);
    swapParticlePositions(10 * indexRowLength + 2, 11 * indexRowLength + 2);
    updateBindings();
    
    EXPECT_TRUE(areParticlesBound(indexRowLength + 4, indexRowLength));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 4, indexRowLength + 4));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 4, 2));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 4, 2 * indexRowLength + 2));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 2, indexRowLength + 6));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 2, 4));
    EXPECT_TRUE(areParticlesBound(indexRowLength + 2, 2 * indexRowLength + 4));

    EXPECT_TRUE(areParticlesBound(11 * indexRowLength + 2, 10 * indexRowLength));
    EXPECT_TRUE(areParticlesBound(11 * indexRowLength + 2, 10 * indexRowLength + 4));
    EXPECT_TRUE(areParticlesBound(11 * indexRowLength + 2, 9 * indexRowLength + 2));
    EXPECT_TRUE(areParticlesBound(11 * indexRowLength + 2, 11 * indexRowLength + 2));
    EXPECT_TRUE(areParticlesBound(10 * indexRowLength + 2, 11 * indexRowLength));
    EXPECT_TRUE(areParticlesBound(10 * indexRowLength + 2, 11 * indexRowLength + 4));
    EXPECT_TRUE(areParticlesBound(10 * indexRowLength + 2, 12 * indexRowLength + 2));
}
