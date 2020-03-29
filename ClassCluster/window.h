#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

class Picture;
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

    public slots:
        void run();
        void restart();
        void move();
        void speedChanged();
        void draw(double *XP, double *YP, double *ZP, int N);

    protected:
        void closeEvent(QCloseEvent *event);

    private:
        void stopCalc();

        QLineEdit *StepE, *EnE, *Speed;
        QPushButton *Start, *Restart, *End, *Rotate, *Move;
        Picture *Pict;
        Calculation *Calc;
};

#endif // WINDOW_H
