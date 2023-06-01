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
        Potential closestTwo, nextTwo, remaining, angular;
        QString dataDir(DATA_DIRECTORY);
        closestTwo.readData(dataDir + "/ClosestTwo.pot");
        nextTwo.readData(dataDir + "/NextTwo.pot");
        remaining.readData(dataDir + "/Remaining.pot");
        angular.readData(dataDir + "/Angular.pot");
        struc[Calculation::ClosestTwo].pot = &closestTwo;
        struc[Calculation::NextTwo].pot = &nextTwo;
        struc[Calculation::Remaining].pot = &remaining;
        struc[Calculation::Angular].pot = &angular;
        Calc = new Calculation(struc);
    }

    void TearDown()
    {

    }

protected:

    bool bindsParticleTo(Particle* P1, Particle* P2) const
    {
        for (int n=0; n < P1->NB; ++n) if (P1->bound[n].p == P2) return true;
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
             && areParticlesBound(Calc->P + index, Calc->P + index - 40)
             && areParticlesBound(Calc->P + index, Calc->P + index - 38)
             && areParticlesBound(Calc->P + index, Calc->P + index + 38)
             && areParticlesBound(Calc->P + index, Calc->P + index + 40));
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
        return !(Calc->isBindingDoubled(particleIndex));
    }
    
    int getParticleNB(const int particleIndex) const
    {
        return Calc->P[particleIndex].NB;
    }
    
    int getParticleMNB(const int particleIndex) const
    {
        return Calc->P[particleIndex].MNB;
    }
    
    static double getRandom()
    {
        return static_cast<double>(rand()) / RAND_MAX;
    }
    
    double createRandomSpeedDistribution()
    {
        double energy = 0.0;
        for (int n=0; n < Calc->N; ++n)
        {
            Calc->P[n].v.setX(0.5 - getRandom());
            Calc->P[n].v.setY(0.5 - getRandom());
            Calc->P[n].v.setZ(0.5 - getRandom());
            Calc->P[n].v *= (1e4 * getRandom());
            energy += Calc->P[n].v.lengthSquared();
        }
        return 0.5 * energy;
    }

    void setParticles(const int N, const Vector* const R, const Vector* const v, const int* const MNB)
    {
        for (int x=0; x < Calc->XS; x++) for (int y=0; y < Calc->YS; y++) for (int z=0; z < Calc->ZS; ++z) Calc->G[x][y][z] = nullptr;
        Calc->N = N;
        Calc->mRandPOF1 = static_cast<double>(N) / (static_cast<double>(RAND_MAX) + 1.0);
        Calc->mRandPOF2 = static_cast<double>(N-1) / (static_cast<double>(RAND_MAX) + 1.0);
        Vector F(double(Calc->XS) / Calc->MaxX, double(Calc->YS) / Calc->MaxY, double(Calc->ZS) / Calc->MaxZ * 1.001);
        for (int n=0; n<N; ++n)
        {
            Calc->initializeParticle(Calc->P[n], 5, 5, R[n], F, true);
            Calc->P[n].v = v[n];
            Calc->P[n].MNB = MNB[n];
        }
    }

    void addParticleBinding(const int n1, const int n2)
    {
        Calc->P[n1].bound[Calc->P[n1].NB].p = Calc->P + n2;
        Calc->P[n2].bound[Calc->P[n2].NB].p = Calc->P + n1;
        Calc->P[n1].bound[Calc->P[n1].NB++].lastDist = Calc->P[n2].bound[Calc->P[n2].NB++].lastDist = (Calc->P[n1].R - Calc->P[n2].R).length();
    }

    Vector getCenterOfMass()
    {
        Vector C;
        for (int n=0; n < Calc->N; ++n) C += Calc->P[n].R;
        return (C / Calc->N);
    }

    Vector getAverageV()
    {
        Vector v;
        for (int n=0; n < Calc->N; ++n) v += Calc->P[n].v;
        return (v / Calc->N);
    }

    Vector getAngularMomentum(const Vector& C)
    {
        Vector L;
        for (int n=0; n < Calc->N; ++n) L += (Calc->P[n].R - C).cross(Calc->P[n].v);
        return L;
    }

    void run(const int maxIteration)
    {
        Calc->mMaxIt = maxIteration;
        Calc->start();
        Calc->wait();
    }

    Calculation* Calc;
};


TEST_F(CalculationTest, CheckParticleBindingInitialisation)
{
    bool even = true;
    int start = 0;
    for (int z=1; z < getPZS() - 1; ++z)
    {
        even = (even ? false : true);
        start += (even ? 38 : 40);
        for (int x=2; x < (even ? getPXS() - 2 : getPXS() - 4); ++x)
        {
            int index = start + x;
            EXPECT_TRUE(areParticleBindingsCorrectlyInitialized(index)) << "index = " << index;
        }
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
    double kinEnergy = createRandomSpeedDistribution();
    Calc->setEnergy(-0.7 * kinEnergy);
    EXPECT_NEAR(0.3 * kinEnergy, Calc->getKineticEnergy(), 1e-5);
}

TEST_F(CalculationTest, Triple)
{
    Vector R[3] = {Vector(32.0, 40.0, 40.0), Vector(40.0, 41.0, 40.0), Vector(48.0, 40.0, 40.0)}, v[3];
    int MNB[3] = {1, 2, 1};
    setParticles(3, R, v, MNB);
    addParticleBinding(0, 1);
    addParticleBinding(1, 2);
    double EB = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CB = getCenterOfMass(), vB = getAverageV(), LB = getAngularMomentum(CB);
    run(1000);
    double EA = Calc->getKineticEnergy() + Calc->getPotentialEnergy();
    Vector CA = getCenterOfMass(), vA = getAverageV(), LA = getAngularMomentum(CA);
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
