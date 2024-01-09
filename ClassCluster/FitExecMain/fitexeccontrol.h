#ifndef FITEXECCONTROL_H
#define FITEXECCONTROL_H


#include <QObject>
#include <QMutex>

#include "potstruct.h"
#include "Calculation.h"


class Calculation;
class CalculationTestHelper;
class Vector;


class FitExecControl : public QObject
{
    Q_OBJECT

public:
    FitExecControl();
    ~FitExecControl();

private slots:
    void printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev);
    void initInstance(int instanceId);
    void calculationStopped(int instanceId);

private:
    void saveResults();
    void addToNullDiff(const double value);
    void addToPairDiff(const double value1, const double value2);
    void adddPdRErrorRatio(const int instanceId, const int iteration, const CalculationTestHelper& helper, const Vector* const a, const double lastOE);

    bool stopped[6] = {false, false, false, false, false, false};
    //double max[6];
    PotStruct struc[Calculation::NumPot];
    double startAngles[93], **angles, **dPdRErrorRatio, **energyDiffs, lastE[6], pairDiffSum = 0.0, maxPairDiff = 0.0, nullDiffSum = 0.0, maxNullDiff = 0.0, *nullPot = nullptr, ecalcFQS = 0.0;
    int currentIndex = -1, instanceIndex[6], maxIteration[6], numPairDiff = 0, numNullDiff = 0, errorCount = 0;
    Calculation* Calc[6];
    QMutex mutex;
};

#endif // FITEXECCONTROL_H
