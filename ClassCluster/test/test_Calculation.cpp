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
#include "potstruct.h"
#include "potential.h"

#include <algorithm>


class CalculationTest : public ::testing::Test
{
public:
    CalculationTest()
        : Calc(nullptr)
    {

    }

    ~CalculationTest()
    {
        if (nullptr != Calc) delete Calc;
    }

    void SetUp()
    {
        PotStruct struc[Calculation::NumPot];
        Potential closestTwo, nextTwo, remaining, secondOrder;
        QString dataDir(DATA_DIRECTORY);
        closestTwo.readData(dataDir + "/ClosestTwo.pot");
        nextTwo.readData(dataDir + "/NextTwo.pot");
        remaining.readData(dataDir + "/Remaining.pot");
        secondOrder.readData(dataDir + "/SecondOrder.pot");
        struc[Calculation::ClosestTwo].pot = &closestTwo;
        struc[Calculation::NextTwo].pot = &nextTwo;
        struc[Calculation::Remaining].pot = &remaining;
        struc[Calculation::SecondOrder].pot = &secondOrder;
        Calc = new Calculation(struc);
    }

    void TearDown()
    {

    }

protected:

    bool bindsParticleTo(Particle* P1, Particle* P2) const
    {
        for (int n=0; n<4; ++n) if (P1->bound[n].p == P2) return true;
        return false;
    }

    bool areParticlesBound(Particle* P1, Particle* P2) const
    {
        return (bindsParticleTo(P1, P2) && bindsParticleTo(P2, P1));
    }

    bool areParticlesBound(const int index1, const int index2) const
    {
        return areParticlesBound(Calc->P + index1, Calc->P + index2);
    }

    bool areParticleBindingsCorrectlyInitialized(int index) const
    {
        return (areParticlesBound(Calc->P + index, Calc->P + index - 2)
             && areParticlesBound(Calc->P + index, Calc->P + index + 2)
             && areParticlesBound(Calc->P + index, Calc->P + index - 2 * Calc->PXS)
             && areParticlesBound(Calc->P + index, Calc->P + index + 2 * Calc->PXS));
    }

    void swapParticlePositions(const int index1, const int index2)
    {
        std::swap(Calc->P[index1].R, Calc->P[index2].R);
    }

    void updateBindings()
    {
        Calc->UpdateBindings();
    }

    int getPXS() const
    {
        return Calc->PXS;
    }

    int getPZS() const
    {
        return Calc->PZS;
    }
    
    double getDist(const int index1, const int index2) const
    {
        return Calculation::dist(Calc->P + index1, Calc->P + index2);
    }
    
    int getBoundParticleIndex(const int particleIndex, const int bindingIndex) const
    {
        return Calc->P[particleIndex].bound[bindingIndex].p - Calc->P;
    }
    
    int getNumParticles() const
    {
        return Calc->N;
    }
    
    bool isNoBindingDoubled(const int particleIndex) const
    {
       for (int i=0; i < Calc->P[particleIndex].NB - 1; ++i) for (int j=i+1; j < Calc->P[particleIndex].NB; ++j)
           if (Calc->P[particleIndex].bound[i].p == Calc->P[particleIndex].bound[j].p) return false;
        return true;
    }
    
    int getParticleNB(const int particleIndex) const
    {
        return Calc->P[particleIndex].NB;
    }
    
    int getParticleMNB(const int particleIndex) const
    {
        return Calc->P[particleIndex].MNB;
    }

    Calculation* Calc;
};

TEST_F(CalculationTest, CheckParticleBindingInitialisation)
{
    for (int z=1; z < getPZS() - 1; ++z) for (int x=1; x < getPXS() - 1; ++x)
    {
        int index = 2 * (getPXS() * z + x);
        EXPECT_TRUE(areParticleBindingsCorrectlyInitialized(index)) << "index = " << index;
        EXPECT_TRUE(areParticleBindingsCorrectlyInitialized(index + 1)) << "index + 1 = " << (index + 1);
    }
    for (int n=0; n < getNumParticles(); ++n) EXPECT_EQ(getParticleNB(n), getParticleMNB(n)) << "index = " << n;
}

TEST_F(CalculationTest, CheckParticleBindingUpdates)
{
    int indexRowLength = 2 * getPXS();
    swapParticlePositions(indexRowLength + 2, indexRowLength + 4);
    swapParticlePositions(10 * indexRowLength + 2, 11 * indexRowLength + 2);
    updateBindings();

    for (int i=0; i < getNumParticles(); ++i)
    {
        EXPECT_TRUE(isNoBindingDoubled(i));
        
        for (int j=0; j < getParticleNB(i); ++j)
        {
            int k = getBoundParticleIndex(i, j);
            if (k<0) break;
            
            EXPECT_TRUE(areParticlesBound(i, k));
            EXPECT_LT(getDist(i, k), 6.0);           
        }
    }
}

TEST_F(CalculationTest, CalcEndpointsOfEnergyDefinitionAxis)
{
    Vector end1, end2, direction(1.0, 1.0, 1.0);
    Calc->CalcEndpointsOfEnergyDefinitionAxis(210, direction, end1, end2);
    EXPECT_EQ(0.0, end1.X());
    EXPECT_NEAR(15.1716, end1.Y(), 1e-4);
    EXPECT_EQ(0.0, end1.Z());
    EXPECT_NEAR(64.8284, end2.X(), 1e-4);
    EXPECT_EQ(80.0, end2.Y());
    EXPECT_NEAR(64.8284, end2.Z(), 1e-4);
}
