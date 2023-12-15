#include "fitexeccontrol.h"

#include <QCoreApplication>
#include <QTextStream>

#include <algorithm>

#include "particle.h"
#include "potential.h"
#include "DiagWindow.h"
#include "CalculationTestHelper.h"
#include "datasortfunctor.h"
#include "utils.h"
#include "heapsort.h"


FitExecControl::FitExecControl()
{
    QString dataDir(DATA_DIRECTORY);
    struc[Calculation::ClosestTwo].InitAsOwner(dataDir + "/ClosestTwo.pot");
    struc[Calculation::NextTwo].InitAsOwner(dataDir + "/NextTwo.pot");
    struc[Calculation::Remaining].InitAsOwner(dataDir + "/Remaining.pot");
    struc[Calculation::Angular].InitAsOwner(dataDir + "/Angular.pot");
    for (int i=0; i < 93; ++i) startAngles[i] = M_PI * (i+8) / 101;
    // angles = Create(93, 1000);
    energyDiffs = Create(93, 1000);
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
            case 0:Destroy
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

FitExecControl::~FitExecControl()
{
    if (nullptr != nullPot) delete[] nullPot;
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
    if (currentIndex >= 93)
    {
        stopped[instanceId] = true;
        return;
    }
    Calc[instanceId] = new Calculation(struc);
    Calc[instanceId]->setInstanceId(instanceId);
    CalculationTestHelper helper(Calc[instanceId]);
    if (nullptr == nullPot)
    {
        const int numPoints = helper.getNumPotentialPoints();
        nullPot = new double[numPoints];
        std::fill(nullPot, nullPot + numPoints, 0.0);
    }
    double currentAngle = 0.5 * startAngles[instanceIndex[instanceId]];
    static const double radius = 8.0, center = 40.0;
    double sinAngR = radius * sin(currentAngle), cosAngR = radius * cos(currentAngle);
    printf("Starting Calc[%d] with index=%d and sinAngR=%f, cosAngR=%f\n", instanceId, instanceIndex[instanceId], sinAngR, cosAngR);
    Vector R[3] = {Vector(center - sinAngR, center + cosAngR, center), Vector(center, center, center), Vector(center + sinAngR, center + cosAngR, center)};
    Vector v[3] = {Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0)};
    int MNB[3] = {1, 2, 1};
    helper.setParticles(3, R, v, MNB);
    helper.addParticleBinding(0, 1);
    helper.addParticleBinding(1, 2);
    lastE[instanceId] = Calc[instanceId]->getPotentialEnergy() + Calc[instanceId]->getKineticEnergy();
    connect(Calc[instanceId], SIGNAL(CalcState(int, int, double, double)), this, SLOT(printCalcState(int, int, double, double)), Qt::DirectConnection);
    connect(Calc[instanceId], SIGNAL(finished()), Calc[instanceId], SLOT(emitStopped()));
    connect(Calc[instanceId], SIGNAL(Stopped(int)), this, SLOT(calculationStopped(int)));
    Calc[instanceId]->start();
}

void FitExecControl::addToNullDiff(const double value)
{
    const double diff = abs(value);
    if (diff > maxNullDiff) maxNullDiff = diff;
    nullDiffSum += diff;
    ++numNullDiff;
}

void FitExecControl::addToPairDiff(const double value1, const double value2)
{
    const double pairDiff = abs(value2 - value1);
    if (pairDiff > maxPairDiff) maxPairDiff = pairDiff;
    pairDiffSum += pairDiff;
    ++numPairDiff;
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
    if (iteration  < maxIt)
    {
        CalculationTestHelper helper(Calc[instanceId]);
        const double currentEnergy = Calc[instanceId]->getPotentialEnergy() + Calc[instanceId]->getKineticEnergy();
        energyDiffs[instanceIndex[instanceId]][iteration] = (currentEnergy - lastE[instanceId]); // / helper.getSpeedSum();
        lastE[instanceId] = currentEnergy;
        // angles[instanceIndex[instanceId]][iteration] = helper.getBindingAngle(0, 1, 2);
        double U(0.0);
        Vector a[3];
        helper.replacePotential(Calculation::Remaining, nullPot, nullPot);
        if (true == helper.getU(0, 2, U, nullptr, 3, a, false))
        {
            addToNullDiff(helper.getPositionDifference(0, 1).dot(a[0]));
            addToNullDiff(helper.getPositionDifference(2, 1).dot(a[2]));
            addToPairDiff(a[2].X(), -1.0 * a[2].Y());
            addToPairDiff(a[0].Y(), a[2].Y());
            addToNullDiff(a[0].Z());
            addToNullDiff(a[1].X());
            addToPairDiff(a[0].Y() + a[2].Y(), -1.0 * a[1].Y());
            addToNullDiff(a[1].Z());
            addToNullDiff(a[2].Z());
        }
        else ++errorCount;
        helper.resetPotential(Calculation::Remaining);
    }
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
    for (int n = maxIteration[instanceId]; n < 1000; ++n)
    {
        energyDiffs[instanceIndex[instanceId]][n] = 0.0;
        // angles[instanceIndex[instanceId]][n] = 0.0;
    }
    initInstance(instanceId);
    for (int i=0; i<6; ++i) if (!stopped[i]) return;
    // for (int i=0; i<6; ++i) printf("instance%d: max=%g\n", i, max[i]);
    saveResults();
    Destroy(energyDiffs, 93);
    // Destroy(angles, 93);
    QCoreApplication::exit(0);
}

void FitExecControl::saveResults()
{
    double FQS = 0.0;
    for (int n=0; n<93; ++n) for (int i=0; i<1000; ++i) FQS += energyDiffs[n][i] * energyDiffs[n][i];
    printf("numError=%d, numPairDiff=%d, numNullDiff=%d, maxPairDiff=%g, avPairDiff=%g, maxNullDiff=%g, avNullDiff=%g\n", errorCount, numPairDiff, numNullDiff, maxPairDiff,
           pairDiffSum / numPairDiff, maxNullDiff, nullDiffSum / numNullDiff);
    printf("FQS=%g\n", FQS);

    /*double **drawData = Create(100000, 2);
    DataSortFunctor sorter(angles, 93, 1000);
    int *sort = utils::heapSort(sorter, 100000);
    int **sortData = sorter.getResult(sort);
    DiagWindow* resultWindow = new DiagWindow;
    for (int n=0; n<100000; ++n)
    {
        drawData[n][0] = angles[sortData[n][0]][sortData[n][1]];
        drawData[n][1] = energyDiffs[sortData[n][0]][sortData[n][1]];
    }
    Destroy(sortData, 100000);
    resultWindow->setData(drawData, 100000);
    resultWindow->show();*/

    /*QFile resultsFile("AngularDeviations.CSV");
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
    }*/
}
