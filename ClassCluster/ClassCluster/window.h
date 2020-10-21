#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

#include "Calculation.h"

class Picture;
class Particle;
class PotStruct;


class Window : public QWidget
{
    Q_OBJECT

    public:
        Window(PotStruct* PotSs = nullptr);
        ~Window();
        void start();
        void stop();
        void reset();
        void move();
        void rotate();
        void triggerSnapShot();
        void restoreSnapShot(bool &isMoving);
        double getEnergy() const;
        double setEnergy(const double E);
        void setPotentialRangeScale(const double newScale);
        void setSpeed(const double newSpeed);
        void setStepSize(const double size);
        void setPotential(const Calculation::PotRole role, PotStruct &Pot);
        void stopCalc();
        bool isRunning() const;
        bool isMoving() const;
        double getRe() const;

    private slots:
        void draw(double *XP, double *YP, double *ZP, int N);
        void writeSnapShot(Particle* P, int N);

    protected:
        void closeEvent(QCloseEvent *event);

    private:

        Picture *Pict;
        Calculation *Calc;
};

#endif // WINDOW_H
