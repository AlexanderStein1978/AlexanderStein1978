#ifndef FITEXECCONTROL_H
#define FITEXECCONTROL_H

#include <QObject>


class Calculation;


class FitExecControl : public QObject
{
    Q_OBJECT

public:
    FitExecControl();

private slots:
    void printCalcState(int instanceId, int iteration, double currentYCenterDev, double maxYCenterDev);

private:
    bool stopped[6] = {false, false, false, false, false, false};
    double max[6];
    Calculation* Calc[6];
};

#endif // FITEXECCONTROL_H
