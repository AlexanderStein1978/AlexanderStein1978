#ifndef FITEXECCONTROL_H
#define FITEXECCONTROL_H


#include <QObject>
#include <QMutex>

#include "potstruct.h"
#include "Calculation.h"


class Calculation;


class FitExecControl : public QObject
{
    Q_OBJECT

public:
    FitExecControl();

private slots:
    void printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev);
    void initInstance(int instanceId);

private:
    void saveResults();

    bool stopped[6] = {false, false, false, false, false, false};
    //double max[6];
    PotStruct struc[Calculation::NumPot];
    double angles[100], **results, startE[6];
    int currentIndex = -1, instanceIndex[6];
    Calculation* Calc[6];
    QMutex mutex;
};

#endif // FITEXECCONTROL_H
