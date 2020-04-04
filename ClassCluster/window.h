#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

class Picture;
class Particle;
class Calculation;

class QLineEdit;
class QCloseEvent;
class QPushButton;


class Window : public QWidget
{
    Q_OBJECT

    public:
        Window();
        ~Window();

    private slots:
        void run();
        void restart();
        void move();
        void speedChanged();
        void draw(double *XP, double *YP, double *ZP, int N);
        void writeSnapShot(Particle* P, int N);

    protected:
        void closeEvent(QCloseEvent *event);

    private:
        void stopCalc();

        QLineEdit *StepE, *EnE, *Speed;
        QPushButton *Start, *Restart, *SnapShot, *Rotate, *Move;
        Picture *Pict;
        Calculation *Calc;
};

#endif // WINDOW_H
