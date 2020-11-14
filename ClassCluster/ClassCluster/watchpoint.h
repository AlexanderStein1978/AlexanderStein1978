#ifndef WATCHPOINT_H
#define WATCHPOINT_H


class WatchStep;


class WatchPoint
{
public:
    WatchPoint(const int iNumSteps, const int iNumParticles);
    ~WatchPoint();
    void reset();
    double get(const int step, const int particleIndex, const int coordinate) const;
    void set(const int step, const int particleIndex, const double x, const double y, const double z);
    void setSum(const int step, const double X, const double Y, const double Z);
    double getSumX(const int step) const;
    double getSumY(const int step) const;
    double getSumZ(const int step) const;

private:
    WatchStep** steps;
    const int numSteps;
};

#endif // WATCHPOINT_H
