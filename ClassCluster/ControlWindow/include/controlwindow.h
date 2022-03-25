#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H


#include <QWidget>
#include <QString>

class QLineEdit;
class QLabel;
class QCloseEvent;
class QPushButton;
class QComboBox;
class QTextStream;

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
    void Init(QString& data);
    void run();
    void restart();
    void move();
    void rotate();
    void speedChanged();
    void writeSnapShot();
    void restoreSnapShot();
    void showParticleWatchWindow();
    void saveSettings();
    void plotClosing();
    void EChanged();
    void TChanged();
    void ValueChanged();
    void EnergyRelevantValueChanged();
    void UpdateEnergies(double kinteticEnergy, double totalEnergy);
    void networkSelectionChanged(int index);
    void connectToServer();
    void setIsRunning(bool isRunning);

    inline void connectionEstablished()
    {
        setConnectionStatus(true);
    }

    void disconnected()
    {
        setConnectionStatus(false);
    }

protected:
    void focusInEvent(QFocusEvent *event) override;

private:
    enum NetworkSelections{NetSelLocal, NetSelServer, NetSelClient};

    void Init(QTextStream& inStream);
    void Serialize(QTextStream& outStream);
    void prepareWindow();
    void start();
    bool stopIfItsRunning();
    void setConnectionStatus(bool connected);

    Window* window;
    QLineEdit *StepE, *EnE, *TEdit, *Speed, *PotRangeScaleEdit, *LayerDistanceEdit, *IpAddressEdit;
    QLabel *KineticEnergyLabel, *PotentialEnergyLabel, *TotalEnergyLabel, *ConnectionStatus;
    QComboBox *NetworkSelection;
    QPushButton *Start, *Restart, *WriteSnapShot, *RestoreSnapShot, *ShowParticleWatchWindow, *Rotate, *Move, *Connect;
    PotControl** PotControls;
    PotentialPlot* Plot;
    MainWindow *MW;
    const QString SettingsFileName, ProgramPath;
};

#endif // CONTROLWINDOW_H
