//
// C++ Implementation: PotFit
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "PotFit.h"
#include "potworker.h"


PotFit::PotFit()
{
    robustWeighting = false;
    UseSVD = true;
    threshhold = 1e-10;
    useLevenbergMarquardt = true;
    MaxIt = -1;
    adjustTEWeightFact = false;
    hCi_cP = 0.0;
    Worker = 0;
}

void PotFit::FitAnaPot(bool nrobustWeighting, bool nUseSVD, double nthreshhold, 
                         bool nuseLevenbergMarquardt, int nMaxIt, bool nadjustTEWeightFact, double nhCi_cP)
{
    robustWeighting = nrobustWeighting;
    UseSVD = nUseSVD;
    threshhold = nthreshhold;
    useLevenbergMarquardt = nuseLevenbergMarquardt;
    MaxIt = nMaxIt;
    adjustTEWeightFact = nadjustTEWeightFact;
    hCi_cP = nhCi_cP;
    FitPot = analyticalPotential;
    restart = false;
    ToDo = PotentialFit;
}

void PotFit::run()
{
    int NumWFPoints = NumPoints;
    setTerminationEnabled(true);
    if (Worker != 0) 
    {
        switch (ToDo)
        {
            case PotentialFit:
                switch (FitPot)
                {
                    case analyticalPotential:
                        Worker->FitAnaPot(NumWFPoints, robustWeighting, UseSVD, threshhold, useLevenbergMarquardt, MaxIt, adjustTEWeightFact, hCi_cP,
                                          restart);
                        break;
                    case SplinePotential:
                        Worker->FitSplinePot(NumWFPoints, MaxIt, threshhold);
                        break;
                    default:
                        printf("Potentials of type %d currently do not be fitted in a separate thread!", FitPot);
                        break;
                }
                break;
            case MonteCarloSimulation:
                Worker->MonteCarloSimIteration(MCST, MCSUncFact);
                break;
            case calcFQS:
                Worker->getFQS(NumWFPoints, false, 0, true, SFQSRad);
                break;
        }
    }
}

void PotFit::setWorker(PotWorker* nWorker)
{
    Worker = nWorker;
}
