#include "fitexeccontrol.h"

#include <QCoreApplication>
#include <QTextStream>

#include "particle.h"
#include "potential.h"
#include "CalculationTestHelper.h"
#include "utils.h"


FitExecControl::FitExecControl()
{
    QString dataDir(DATA_DIRECTORY);
    struc[Calculation::ClosestTwo].InitAsOwner(dataDir + "/ClosestTwo.pot");
    struc[Calculation::NextTwo].InitAsOwner(dataDir + "/NextTwo.pot");
    struc[Calculation::Remaining].InitAsOwner(dataDir + "/Remaining.pot");
    struc[Calculation::Angular].InitAsOwner(dataDir + "/Angular.pot");
    for (int i=0; i<100; ++i) angles[i] = M_PI * (i+1) / 101;
    results = Create(100, 1000);
    for (int i=0; i<6; ++i)
    {
        /*int i = 0;
        currentIndex = 98;*/

        Calc[i] = nullptr;
        /*Calc[i]->setLayerDistance(20.0);
        Calc[i]->setEnergy(200000.0);
        Calc[i]->move();
        */
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
        /*switch(i)
        {
            case 0:
                Calc[i]->setStepSize(5e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE1_200KE.csv");
                break;
            case 1:
                Calc[i]->setStepSize(2e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE2_200KE.csv");
                break;
            case 2:
                Calc[i]->setStepSize(1e-4);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE3_200KE.csv");
                break;
            case 3:
                Calc[i]->setStepSize(5e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE4_200KE.csv");
                break;
            case 4:
                Calc[i]->setStepSize(2e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE5.csv");
                break;
            case 5:
                Calc[i]->setStepSize(1e-5);
                Calc[i]->setEnergyCsvLogFilename("CSVEnergyLogFE6.csv");
                break;
        }*/
        initInstance(i);
    }
}

void FitExecControl::initInstance(int instanceId)
{
    maxIteration[instanceId] = -1;
    if (nullptr != Calc[instanceId])
    {
        disconnect(Calc[instanceId], SIGNAL(CalcState(int, int, double, double)), this, SLOT(printCalcState(int, int, double, double)));
        disconnect(Calc[instanceId], SIGNAL(Stopped(int)), this, SLOT(calculationStopped(int)));
        delete Calc[instanceId];
    }
    mutex.lock();
    instanceIndex[instanceId] = ++currentIndex;
    mutex.unlock();
    if (currentIndex >= 100)
    {
        stopped[instanceId] = true;
        return;
    }
    Calc[instanceId] = new Calculation(struc);
    Calc[instanceId]->setInstanceId(instanceId);
    CalculationTestHelper helper(Calc[instanceId]);
    double currentAngle = 0.5 * angles[instanceIndex[instanceId]];
    static const double radius = 8.0, center = 40.0;
    double sinAngR = radius * sin(currentAngle), cosAngR = radius * cos(currentAngle);
    printf("Starting Calc[%d] with index=%d and sinAngR=%f, cosAngR=%f\n", instanceId, instanceIndex[instanceId], sinAngR, cosAngR);
    Vector R[3] = {Vector(center - sinAngR, center + cosAngR, center), Vector(center, center, center), Vector(center + sinAngR, center + cosAngR, center)};
    Vector v[3] = {Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0)};
    int MNB[3] = {1, 2, 1};
    helper.setParticles(3, R, v, MNB);
    helper.addParticleBinding(0, 1);
    helper.addParticleBinding(1, 2);
    startE[instanceId] = Calc[instanceId]->getPotentialEnergy() + Calc[instanceId]->getKineticEnergy();
    connect(Calc[instanceId], SIGNAL(CalcState(int, int, double, double)), this, SLOT(printCalcState(int, int, double, double)), Qt::DirectConnection);
    connect(Calc[instanceId], SIGNAL(finished()), Calc[instanceId], SLOT(emitStopped()));
    connect(Calc[instanceId], SIGNAL(Stopped(int)), this, SLOT(calculationStopped(int)));
    Calc[instanceId]->start();
}


void FitExecControl::printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev)
{
    //printf("instance%d: it=%d, currentDev=%g, maxDev=%g\n", instanceId, iteration, currentYCenterDev, maxYCenterDev);
    static const int maxIt = 1000;
    /*switch(instanceId)
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
    }*/
    maxIteration[instanceId] = iteration;
    if (iteration  < maxIt) results[instanceIndex[instanceId]][iteration] = Calc[instanceId]->getPotentialEnergy() + Calc[instanceId]->getKineticEnergy() - startE[instanceId];
    else
    {
        //max[instanceId] = maxYCenterDev;
        if (!stopped[instanceId])
        {
            Calc[instanceId]->stop();
            // Calc[instanceId]->wait();
            // printf("Calculation with instanceIndex[%d]=%d finished.", instanceId, instanceIndex[instanceId]);
            // initInstance(instanceId);
        }

    }
}

void FitExecControl::calculationStopped(int instanceId)
{
    for (int n = maxIteration[instanceId]; n < 1000; ++n) results[instanceIndex[instanceId]][n] = 0.0;
    initInstance(instanceId);
    for (int i=0; i<6; ++i) if (!stopped[i]) return;
    // for (int i=0; i<6; ++i) printf("instance%d: max=%g\n", i, max[i]);
    saveResults();
    Destroy(results, 100);
    QCoreApplication::exit(0);
}

void FitExecControl::saveResults()
{
    QFile resultsFile("AngularDeviations.CSV");
    resultsFile.open(QIODevice::WriteOnly);
    QTextStream S(&resultsFile);
    double AF = 180.0 / M_PI;
    for (int i=0; i<100; ++i) S << '\t' << QString::number(AF * angles[i], 'f', 2);
    S << '\n';
    for (int i=0; i < 1000; ++i)
    {
        S << QString::number(i);
        for (int j=0; j < 100; ++j) S << '\t' << QString::number(results[j][i]);
        S << '\n';
    }
}
