#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

class Picture;
class Particle;
class Calculation;


class Window : public QWidget
{
    Q_OBJECT

    public:
        Window();
        ~Window();
        void start();
        void stop();
        void reset();
        void move();
        void rotate();
        void triggerSnapShot();
        void restoreSnapShot();
        double getEnergy() const;
        double setEnergy(const double E);
        void setSpeed(const double newSpeed);
        void setStepSize(const double size);
        void stopCalc();
        bool isRunning() const;
        bool isMoving() const;

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
