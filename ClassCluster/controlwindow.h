#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H


#include <QWidget>


class QLineEdit;
class QCloseEvent;
class QPushButton;

class Window;


class ControlWindow : public QWidget
{
    Q_OBJECT

public:
    ControlWindow();
    ~ControlWindow();

private slots:
    void run();
    void restart();
    void move();
    void rotate();
    void speedChanged();
    void writeSnapShot();
    void restoreSnapShot();

protected:
    void closeEvent(QCloseEvent *event);

private:
    Window* window;

    QLineEdit *StepE, *EnE, *Speed;
    QPushButton *Start, *Restart, *WriteSnapShot, *RestoreSnapShot, *Rotate, *Move;

};

#endif // CONTROLWINDOW_H
