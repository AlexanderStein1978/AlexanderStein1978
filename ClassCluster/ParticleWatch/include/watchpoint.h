//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef WATCHPOINT_H
#define WATCHPOINT_H


class WatchStep;
class Vector;


class WatchPoint
{
public:
    WatchPoint(const int iNumSteps, const int iNumParticles);
    ~WatchPoint();
    void reset();
    double get(const int step, const int particleIndex, const int coordinate) const;
    void set(const int step, const int particleIndex, const Vector& r);
    void setSum(const int step, const Vector& R);
    double getSumX(const int step) const;
    double getSumY(const int step) const;
    double getSumZ(const int step) const;

private:
    WatchStep** steps;
    const int numSteps;
};

#endif // WATCHPOINT_H
