#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

#include "Calculation.h"

class Picture;
class Particle;
class PotStruct;
class WatchPoint;


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
        double getPotentialEnergy() const;
        double getKineticEnergy() const;
        double setKineticEnergy(const double newT);
        void setPotentialRangeScale(const double newScale);
        void setSpeed(const double newSpeed);
        void setStepSize(const double size);
        void setPotential(const Calculation::PotRole role, PotStruct &Pot);
        void stopCalc();
        bool isRunning() const;
        bool isMoving() const;
        double getRe() const;
        void setParticleWatchPoint(WatchPoint* point);
        void setParticleWatch(const int indexToWatch);
        int getXDim() const;
        int getNumParticles() const;
        static int getNumSteps();
        void setLayerDistance(const double newDistance);

    signals:
        void EnergiesChanged(double kineticEnergy, double totalEnergy);

    private slots:
        void draw(Vector* Pos, int N);
        void writeSnapShot(Particle* P, int N);

    protected:
        void closeEvent(QCloseEvent *event) override;
        void paintEvent(QPaintEvent *e) override;

    private:
        void destroyData();

        Picture *Pict;
        Calculation *Calc;

        Vector *mPos;
        int mN;
};

#endif // WINDOW_H
