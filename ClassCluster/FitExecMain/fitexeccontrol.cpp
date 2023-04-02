#include "fitexeccontrol.h"

#include <QCoreApplication>

#include "Calculation.h"
#include "particle.h"
#include "potstruct.h"
#include "potential.h"


FitExecControl::FitExecControl()
{
    PotStruct struc[Calculation::NumPot];
    Potential closestTwo, nextTwo, remaining, secondOrder, angular;
    QString dataDir(DATA_DIRECTORY);
    closestTwo.readData(dataDir + "/ClosestTwo.pot");
    nextTwo.readData(dataDir + "/NextTwo.pot");
    remaining.readData(dataDir + "/Remaining.pot");
    secondOrder.readData(dataDir + "/SecondOrder.pot");
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
        switch(i)
        {
            case 0:
                Calc[i]->setWaveStep(10);
            case 1:
                Calc[i]->setWaveStep(12.5);
            case 2:
                Calc[i]->setWaveStep(15);
            case 3:
                Calc[i]->setWaveStep(17.5);
            case 4:
                Calc[i]->setWaveStep(20);
            case 5:
                Calc[i]->setWaveStep(50);
        }
        connect(Calc[i], SIGNAL(CalcState(int, int, double, double)), this, SLOT(printCalcState(int, int, double, double)));
        Calc[i]->start();
    }
}

void FitExecControl::printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev)
{
    printf("instance%d: it=%d, currentDev=%g, maxDev=%g\n", instanceId, iteration, currentYCenterDev, maxYCenterDev);
    if (iteration >= 6000)
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
