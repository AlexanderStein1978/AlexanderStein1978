#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H


#include <QWidget>
#include <QString>

class QLineEdit;
class QCloseEvent;
class QPushButton;

class Window;
class PotControl;
class PotentialPlot;


class ControlWindow : public QWidget
{
    Q_OBJECT

public:
    ControlWindow();
    ~ControlWindow();
    double getRe() const;
    
    inline const QString& getProgramPath() const
    {
        return ProgramPath;
    }

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

    QLineEdit *StepE, *EnE, *Speed, *PotRangeScaleEdit;
    QPushButton *Start, *Restart, *WriteSnapShot, *RestoreSnapShot, *Rotate, *Move;
    PotControl** PotControls;
    PotentialPlot *const Plot;
    const QString SettingsFileName, ProgramPath;
};

#endif // CONTROLWINDOW_H
