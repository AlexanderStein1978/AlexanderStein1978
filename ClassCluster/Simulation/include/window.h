#ifndef WINDOW_H
#define WINDOW_H


#include <QWidget>

#include "Calculation.h"


class QTcpServer;
class QTcpSocket;

class Picture;
class Particle;
class WatchPoint;
class NetworkClient;
class NetworkServer;

struct PotStruct;


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
        void restoreSnapShot(bool &isMoving, const QString& filename);
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
        bool listenAsCalculationServer(const QString IpAddress);
        void connectToCalculationServer(const QString IpAddress);
        void switchBackToLocalCalulations();
        void stopListeningAsCalculationServer();
        void sendReloadPotentials();
        void writeSnapShot(QString fileName);
        void emitReloadPotentials();
        void copyDataIfNew(QByteArray& data, const QByteArray& sendCommand);
        char* getDrawingDataToFill(const int N);
        void flagsReceived(char flags);
        void SendSettings(const QByteArray& data);
        void SendPotential(const QByteArray& data);
        void SendGetSettingsRequest();
        void SettingsReceived(const QByteArray& data);
        void PotentialReceived(const QByteArray& data);
        void GetSettingsRequestReceived();

    signals:
        void EnergiesChanged(double kineticEnergy, double totalEnergy);
        void ReloadPotentials();
        void IsConnectedToServer();
        void ConnectionFailed();
        void IsRunning(bool running);
        void SendSettings();
        void ReceivedSetting(const QByteArray& data);
        void ReceivedPotential(const QByteArray& data);

    public slots:
        void updateRemoteEnergies(const double kineticEnergy, const double totalEnergy);

    private slots:
        void draw(Vector* Pos, int N);
        void writeSnapShot(Particle* P, int N);
        void newClientConnection();

    protected:
        void closeEvent(QCloseEvent *event) override;
        void paintEvent(QPaintEvent *e) override;

    private:
        void destroyData();
        void writeSnapShot(Particle* P, int N, QString& fileName);

        Picture *Pict;
        Calculation *Calc;

        QTcpServer* mServer;
        NetworkClient* mNetworkClient;
        NetworkServer* mNetworkServer;
        bool mIsRemoteRunning, mDataIsNew, mWaitingForData;
        QMutex mDataMutex;

        Vector *mPos;
        int mN, mNumRemoteParticles;
        double mRemoteKineticEnergy, mRemoteTotalEnergy;
        QString mRemoteSnapShotFileName;
        Particle* mRemoteSnapShotParticles;
};

#endif // WINDOW_H
