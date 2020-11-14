#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H


#include <QWidget>
#include <QString>

class QLineEdit;
class QCloseEvent;
class QPushButton;

class Window;
class PotControl;
class Potential;
class PotentialPlot;
class MainWindow;


class ControlWindow : public QWidget
{
    Q_OBJECT

public:
    ControlWindow(MainWindow *const mw);
    ~ControlWindow();
    double getRe() const;
    void showPotential(Potential* const pot, const bool plot);
    
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
    void showParticleWatchWindow();
    void saveSettings();

protected:
    void focusInEvent(QFocusEvent *event) override;

private:
    void prepareWindow();

    Window* window;
    QLineEdit *StepE, *EnE, *Speed, *PotRangeScaleEdit;
    QPushButton *Start, *Restart, *WriteSnapShot, *RestoreSnapShot, *ShowParticleWatchWindow, *Rotate, *Move;
    PotControl** PotControls;
    PotentialPlot* Plot;
    MainWindow *MW;
    const QString SettingsFileName, ProgramPath;
};

#endif // CONTROLWINDOW_H
