#include "fitexeccontrol.h"

#include <QCoreApplication>

#include "Calculation.h"
#include "particle.h"
#include "potstruct.h"
#include "potential.h"


FitExecControl::FitExecControl()
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
    for (int i=0; i<6; ++i)
    {
        Calc[i] = new Calculation(struc);
        Calc[i]->setLayerDistance(20.0);
        Calc[i]->setEnergy(400000.0);
        Calc[i]->move();
        Calc[i]->setInstanceId(i);
        /*switch(i)
        {
            case 0:
                Calc[i]->setWaveStep(10);
                break;
            case 1:
                Calc[i]->setWaveStep(12.5);
                break;
            case 2:
                Calc[i]->setWaveStep(15);
                break;
            case 3:
                Calc[i]->setWaveStep(17.5);
                break;
            case 4:
                Calc[i]->setWaveStep(20);
                break;
            case 5:
                Calc[i]->setWaveStep(50);
                break;
        }*/
        switch(i)
        {
            case 0:
                Calc[i]->setStepSize(5e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE1.csv");
                break;
            case 1:
                Calc[i]->setStepSize(2e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE2.csv");
                break;
            case 2:
                Calc[i]->setStepSize(1e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE3.csv");
                break;
            case 3:
                Calc[i]->setStepSize(5e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE4.csv");
                break;
            case 4:
                Calc[i]->setStepSize(2e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE5.csv");
                break;
            case 5:
                Calc[i]->setStepSize(1e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE6.csv");
                break;
        }

        connect(Calc[i], SIGNAL(CalcState(int, int, double, double)), this, SLOT(printCalcState(int, int, double, double)));
        Calc[i]->start();
    }
}

void FitExecControl::printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev)
{
    printf("instance%d: it=%d, currentDev=%g, maxDev=%g\n", instanceId, iteration, currentYCenterDev, maxYCenterDev);
    int maxIt;
    switch(instanceId)
    {
        case 0:
            maxIt = 2000;
            break;
        case 1:
            maxIt = 5000;
            break;
        case 2:
            maxIt = 10000;
            break;
        case 3:
            maxIt = 20000;
            break;
        case 4:
            maxIt = 50000;
            break;
        case 5:
            maxIt = 100000;
            break;
    }
    if (iteration >= maxIt)
    {
        max[instanceId] = maxYCenterDev;
        if (!stopped[instanceId])
        {
            Calc[instanceId]->stop();
            stopped[instanceId] = true;
        }
        for (int i=0; i<6; ++i) if (!stopped[i]) return;
        for (int i=0; i<6; ++i) printf("instance%d: max=%g\n", i, max[i]);
        QCoreApplication::exit(0);
    }
}
