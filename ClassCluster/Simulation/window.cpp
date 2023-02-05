#include "window.h"
#include "Calculation.h"
#include "Picture.h"
#include "particle.h"
#include "networkclient.h"
#include "networkserver.h"

#include <QCloseEvent>
#include <QPainter>
#include <QRect>
#include <QFileDialog>
#include <QTextStream>
#include <QBoxLayout>
#include <QMessageBox>
#include <QTcpServer>
#include <QTcpSocket>


Window::Window(PotStruct *PotSs) : mServer(nullptr), mNetworkClient(nullptr), mNetworkServer(nullptr), mIsRemoteRunning(false), mDataIsNew(false), mWaitingForData(false), mPos(nullptr),
    mN(0), mNumRemoteParticles(0), mMarkedParticle(-1), mMarkedIndex(-1), mRemoteKineticEnergy(0.0), mRemoteTotalEnergy(0.0), mRemoteSnapShotParticles(nullptr)
{
    Calc = new Calculation(PotSs);
    QBoxLayout *L = new QBoxLayout(QBoxLayout::LeftToRight, this);
    L->addWidget(Pict = new Picture(this));
    int XSize, YSize;
    Calc->getSize(XSize, YSize);
    setMinimumSize(XSize, YSize);
    setMaximumSize(XSize, YSize);

    connect(Calc, SIGNAL(WriteSnapShot(Particle*, int)), this, SLOT(writeSnapShot(Particle*, int)));
    connect(Calc, SIGNAL(PictureChanged(Vector*, int)), this, SLOT(draw(Vector*, int)));
    connect(Calc, SIGNAL(EnergiesChanged(double,double)), this, SLOT(updateRemoteEnergies(double,double)));
}

Window::~Window()
{
    stopCalc();
    delete Calc;
    destroyData();
    if (nullptr != mNetworkClient) delete mNetworkClient;
    if (nullptr != mNetworkServer) delete mNetworkServer;
    if (nullptr != mServer) delete mServer;
    if (nullptr != mRemoteSnapShotParticles) delete[] mRemoteSnapShotParticles;
}

void Window::destroyData()
{
    if (0 != mN) delete[] mPos;
}

bool Window::listenAsCalculationServer(const QString IpAddress)
{
    if (nullptr == mServer) mServer = new QTcpServer;
    if (nullptr != mServer)
    {
        bool result = mServer->listen(QHostAddress(IpAddress), 50000);
        if (result) connect(mServer, SIGNAL(newConnection()), this, SLOT(newClientConnection()));
        return result;
    }
    return false;
}

void Window::stopListeningAsCalculationServer()
{
    if (nullptr != mServer)
    {
        if (nullptr != mNetworkServer)
        {
            delete mNetworkServer;
            mNetworkServer = nullptr;
        }
        disconnect(mServer, SIGNAL(newConnection()), this, SLOT(newClientConnection()));
        delete mServer;
        mServer = nullptr;
    }
}

void Window::newClientConnection()
{
    QTcpSocket* newConnection;
    while ((newConnection = mServer->nextPendingConnection()) != nullptr)
    {
        if (nullptr == mNetworkServer) mNetworkServer = new NetworkServer(this, newConnection);
        else mNetworkServer->NewConnection(newConnection);
    }
    char flags(0x00);
    if (Calc->isRunning())
    {
        mDataIsNew = true;
        mNetworkServer->SendData();
        flags = 0x80;
    }
    mNetworkServer->SendFlags(flags);
}

void Window::connectToCalculationServer(const QString IpAddress)
{
    if (nullptr == mNetworkClient)
    {
        mNetworkClient = new NetworkClient(this);
        connect(mNetworkClient, SIGNAL(ConnectionFailed()), this, SIGNAL(ConnectionFailed()));
        connect(mNetworkClient, SIGNAL(IsConnected()), this, SIGNAL(IsConnectedToServer()));
    }
    mNetworkClient->ConnectToServer(IpAddress);
}

void Window::switchBackToLocalCalulations()
{
    if (nullptr != mNetworkClient)
    {
        delete mNetworkClient;
        mNetworkClient = nullptr;
    }
}

void Window::stopCalc()
{
    if (nullptr != mNetworkClient)
    {
        mNetworkClient->SendCommand(Network::STOP_CALC);
        mIsRemoteRunning = false;
    }
    else if (Calc->isRunning())
    {
        Calc->stop();
        Calc->wait();
        if (nullptr != mNetworkServer) mNetworkServer->SendFlags(0x00);
    }
    emit IsRunning(false);
}

void Window::SendSettings(const QByteArray &data)
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::SET_SETTINGS, data);
}

void Window::SendPotential(const QByteArray &data)
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::SET_POTENTIAL, data);
}

void Window::SendGetSettingsRequest()
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::GET_SETTINGS_AND_POTENTIALS);
}

void Window::SettingsReceived(const QByteArray &data)
{
    emit ReceivedSetting(data);
}

void Window::PotentialReceived(const QByteArray &data)
{
    emit ReceivedPotential(data);
}

void Window::GetSettingsRequestReceived()
{
    emit SendSettings();
}

void Window::flagsReceived(char flags)
{
    emit IsRunning((flags & 0x80) != 0x00);
}


void Window::closeEvent(QCloseEvent* event)
{
    if (nullptr == mNetworkClient) stopCalc();
    event->accept();
}

void Window::draw(Vector *Pos, int N)
{
    mDataMutex.lock();
    Calc->mutex.lock();
    if (mN != N)
    {
        destroyData();
        mN = N;
        mPos = new Vector[N];
    }
    memcpy(mPos, Pos, N * sizeof(Vector));
    mDataIsNew = true;
    Calc->mutex.unlock();
    mDataMutex.unlock();
    if (mMarkedIndex >= 0) mMarkedParticle = Calc->TranslateParticleIndex(mMarkedIndex);
    update();
    if (mWaitingForData && nullptr != mNetworkServer) mNetworkServer->SendData();
}

void Window::copyDataIfNew(QByteArray &data, const QByteArray &sendCommand)
{
    QMutexLocker lock(&mDataMutex);
    if (!mDataIsNew)
    {
        mWaitingForData = true;
        return;
    }
    mWaitingForData = mDataIsNew = false;
    const size_t headerSize = 28;
    const size_t dataSize = mN * 24 + headerSize;
    data.reserve(dataSize);
    data += sendCommand;
    data.resize(dataSize);
    char* dataPtr = data.data();
    quint32 remoteN = mN;
    memcpy(dataPtr + 8, &mRemoteKineticEnergy, 8);
    memcpy(dataPtr + 16, &mRemoteTotalEnergy, 8);
    memcpy(dataPtr + 24, &remoteN, 4);
    memcpy(dataPtr + headerSize, mPos, dataSize - headerSize);
}

char* Window::getDrawingDataToFill(const int N)
{
    if (mN != N)
    {
        destroyData();
        mN = N;
        mPos = new Vector[N];
    }
    return reinterpret_cast<char*>(mPos);
}

void Window::updateRemoteEnergies(const double kineticEnergy, const double totalEnergy)
{
    mRemoteKineticEnergy = kineticEnergy;
    mRemoteTotalEnergy = totalEnergy;
    emit EnergiesChanged(kineticEnergy, totalEnergy);
}

void Window::paintEvent(QPaintEvent *e)
{
    if (0 != mN)
    {
        int w = Pict->width(), he = Pict->height(), col, n, z, MaxX, MaxY, MaxZ, markedX = -10, markedY = -10;
        double ScF;
        Calc->getScales(ScF, MaxX, MaxY, MaxZ);
        double zSc = 255.0 / MaxZ, xc, yc;
        QPainter Paint(Pict->getPixmap());
        Paint.eraseRect(0, 0, w, he);
        for (n=0, z=-100; n < mN; n++)
        {
            if (mPos[n].X() <= 0 || mPos[n].X() >= MaxX || mPos[n].Y() <= 0 || mPos[n].Y() >= MaxY) continue;
            if (int(mPos[n].Z() * zSc) != z || n-1 == mMarkedParticle)
            {
                col = 254 - (z = int(mPos[n].Z() * zSc));
                if (col < 0) col = 0;
                if (col > 254) col = 254;
                Paint.setPen(QColor(col, col, col));
                Paint.setBrush(QColor(col, col, col));
            }
            xc = 5 + mPos[n].X() * ScF;
            yc = 5 + mPos[n].Y() * ScF;
            if (n == mMarkedParticle)
            {
                Paint.setPen(QColor(255, 0, 0));
                Paint.setBrush(QColor(255, 0, 0));
                markedX = xc;
                markedY = yc;
            }
            if (abs(markedX - xc) > 2 || abs(markedY - yc) > 2 || n == mMarkedParticle) Paint.drawEllipse(QRectF(xc - 5.0, yc - 5.0, 10.0, 10.0));
        }
        if (mMarkedParticle >= 0)
        {
            Paint.setPen(QColor(0, 0, 0));
            Vector start, stop;
            if (Calc->IsRotated())
            {
                start = Vector(5 + mAxisStart.X() * ScF, 5 + mAxisStart.Z() * ScF, 0.0);
                stop = Vector(5 + mAxisEnd.X() * ScF, 5 + mAxisEnd.Z() * ScF, 0.0);
            }
            else
            {
                start = Vector(5 + mAxisStart.X() * ScF, 5 + mAxisStart.Y() * ScF, 0.0);
                stop = Vector(5 + mAxisEnd.X() * ScF, 5 + mAxisEnd.Y() * ScF, 0.0);
            }
            Paint.drawLine(start.X(), start.Y(), stop.X(), stop.Y());
            Vector step(0.1 * (stop - start));
            double dl = 10.0 / step.length();
            Vector ortho(dl * step.Y(), -dl * step.X(), 0.0), point(start + step);
            for (int n=1; n<=9; ++n, point += step) Paint.drawLine(point.X() - ortho.X(), point.Y() - ortho.Y(), point.X() + ortho.X(), point.Y() + ortho.Y());
        }
        Pict->update();
    }
    e->accept();
}

int Window::getNumParticles() const
{
    return Calc->getNumParticles();
}

double Window::getRe() const
{
    return Calc->getRe();
}

int Window::getNumSteps()
{
    return Calculation::getNumACalcsPerIt();
}

int Window::getXDim() const
{
    return Calc->getNumXDimParticles();
}

void Window::start()
{
    if (nullptr != mNetworkClient)
    {
        mNetworkClient->SendCommand(mNetworkClient->START);
        mIsRemoteRunning = true;
    }
    else
    {
        Calc->start();
        if (nullptr != mNetworkServer) mNetworkServer->SendFlags(0x80);
    }
    emit IsRunning(true);
}

void Window::StartCommandReceived()
{
    emit ReceivedStartCommand();
}

void Window::stop()
{
    if (nullptr != mNetworkClient)
    {
        mNetworkClient->SendCommand(mNetworkClient->STOP);
        mIsRemoteRunning = false;
    }
    else
    {
        Calc->stop();
        if (nullptr != mNetworkServer) mNetworkServer->SendFlags(0x00);
    }
    emit IsRunning(false);
}

bool Window::isRunning() const
{
    return (mNetworkClient == nullptr ? Calc->isRunning() : mIsRemoteRunning);
}

void Window::reset()
{
    if (nullptr != mNetworkClient)
    {
        mNetworkClient->SendCommand(Network::RESET);
        mIsRemoteRunning = false;
    }
    else
    {
        stopCalc();
        Calc->initialize();
    }
}

void Window::rotate()
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::ROTATE);
    else Calc->rotate();
    mMarkedParticle = Calc->TranslateParticleIndex(mMarkedIndex);
}

void Window::move()
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::MOVE);
    Calc->move();
}

bool Window::isMoving() const
{
    return Calc->getMove();
}

void Window::setSpeed(const double newSpeed)
{
    if (nullptr == mNetworkClient) Calc->setSpeed(newSpeed);
}

double Window::getPotentialEnergy() const
{
    return (nullptr != mNetworkClient ? mRemoteTotalEnergy - mRemoteKineticEnergy : Calc->getPotentialEnergy());
}

double Window::getKineticEnergy() const
{
    return (nullptr != mNetworkClient ? mRemoteKineticEnergy : Calc->getKineticEnergy());
}

double Window::setKineticEnergy(const double T)
{
    if (nullptr != mNetworkClient) return mRemoteKineticEnergy;
    return Calc->setKineticEnergy(T);
}

void Window::setParticleWatchPoint(WatchPoint *point)
{
    if (nullptr != mNetworkClient) QMessageBox::information(this, "ClassCluster", "The feature \"setParticleWatchPoint\" is currently not implemented for remote sessions!");
    else Calc->setParticleWatchPoint(point);
}

void Window::setParticleWatch(const int indexToWatch)
{
    if (nullptr != mNetworkClient)
    {
        QMessageBox::information(this, "ClassCluster", "The feature \"setParticleWatch\" is currently not implemented for remote sessions!");
        return;
    }
    bool CalcWasRunning = Calc->isRunning();
    if (CalcWasRunning)
    {
        Calc->stop();
        Calc->wait();
    }
    Calc->setParticleWatchIndex(indexToWatch);
    Calc->start();
    Calc->wait();
    Calc->setParticleWatchIndex(-1);
    if (CalcWasRunning) Calc->start();
}

void Window::setPotentialRangeScale(const double newScale)
{
    if (nullptr == mNetworkClient) Calc->setPotRangeScale(newScale);
}

void Window::setStepSize(const double size)
{
    if (nullptr == mNetworkClient) Calc->setStepSize(size);
}

void Window::triggerSnapShot()
{
    if (nullptr != mNetworkClient)
    {
        mNetworkClient->SendCommand(Network::TRIGGER_SNAP_SHOT);
        QString FN = QFileDialog::getSaveFileName(this, "Select filename and path for snapshot");
        if (nullptr != FN) mNetworkClient->SendCommand(mNetworkClient->WRITE_SNAP_SHOT, FN);
    }
    Calc->triggerSnapShot();
}

void Window::setPotential(const Calculation::PotRole role, PotStruct &Pot)
{
    if (nullptr == mNetworkClient) Calc->setPotential(role, Pot);
}

void Window::sendReloadPotentials()
{
    if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::RELOAD_POTENTIALS);
}

void Window::emitReloadPotentials()
{
    emit ReloadPotentials();
}

void Window::writeSnapShot(QString fileName)
{
    if (nullptr != mRemoteSnapShotParticles)
    {
        writeSnapShot(mRemoteSnapShotParticles, mNumRemoteParticles, fileName);
        delete[] mRemoteSnapShotParticles;
        mRemoteSnapShotParticles = nullptr;
    }
    else mRemoteSnapShotFileName = fileName;
}

void Window::writeSnapShot(Particle *P, int N)
{
    if (nullptr != mNetworkServer)
    {
        if (mRemoteSnapShotFileName.size() > 0)
        {
            writeSnapShot(P, N, mRemoteSnapShotFileName);
            mRemoteSnapShotFileName.clear();
        }
        else
        {
            if (nullptr == mRemoteSnapShotParticles) mRemoteSnapShotParticles = new Particle[N];
            memcpy(mRemoteSnapShotParticles, P, sizeof(Particle) * N);
            mNumRemoteParticles = N;
        }
    }
    else
    {
        QString FN = QFileDialog::getSaveFileName(this, "Select filename and path for snapshot");
        if (!FN.isEmpty()) writeSnapShot(P, N, FN);
    }
}

void Window::writeSnapShot(Particle *P, int N, QString &fileName)
{
    QFile file(fileName);
    file.open(QIODevice::WriteOnly);
    QTextStream S(&file);
    S << "Speed:\t" << QString::number(Calc->getSpeed(), 'f', 12) << "\n";
    S << "StepSize:\t" << QString::number(Calc->getStepSize(), 'f', 12) << "\n";
    S << (Calc->getMove() ? "is moving\n" : "is not moving\n");
    S << " xp \t yp \t zp \t X \t Y \t Z \t vX \t vY \t vZ \t aaX \t aaY \t aaZ \t lX \t lY \t lZ \t lvX \t lvY \t lvZ \n";
    for (int n=0; n<N; ++n)
        S << P[n].xp << "\t" << P[n].yp << "\t" << P[n].zp << "\t" << QString::number(P[n].R.X(), 'f', 12) << "\t" << QString::number(P[n].R.Y(), 'f', 12) << "\t"
          << QString::number(P[n].R.Z(), 'f', 12) << "\t" << QString::number(P[n].v.X(), 'f', 12) << "\t" << QString::number(P[n].v.Y(), 'f', 12) << "\t"
          << QString::number(P[n].v.Z(), 'f', 12) << "\t" << QString::number(P[n].aa.X(), 'f', 12) << "\t" << QString::number(P[n].aa.Y(), 'f', 12) << "\t"
          << QString::number(P[n].aa.Z(), 'f', 12) << "\t" << QString::number(P[n].lR.X(), 'f', 12) << "\t" << QString::number(P[n].lR.Y(), 'f', 12) << "\t"
          << QString::number(P[n].lR.Z(), 'f', 12) << "\t" << QString::number(P[n].lv.X(), 'f', 12) << "\t" << QString::number(P[n].lv.Y(), 'f', 12) << "\t"
          << QString::number(P[n].lv.Z(), 'f', 12) << "\n";
}

void Window::restoreSnapShot(bool &isMoving)
{
    QString FileN = QFileDialog::getOpenFileName(this, "Select the snapshot to restore");
    if (!FileN.isEmpty())
    {
        if (nullptr != mNetworkClient) mNetworkClient->SendCommand(Network::RESTORE_SNAP_SHOT);
        else restoreSnapShot(isMoving, FileN);
    }
}

void Window::restoreSnapShot(bool &isMoving, const QString& FileN)
{
    int n, N;
    double speed, stepSize;
    bool error = false;
    Calc->initialize();
    QFile file(FileN);
    file.open(QIODevice::ReadOnly);
    QTextStream S(&file);
    QString Buffer = S.readLine();
    QStringList L = Buffer.split('\t');
    if (L.size() != 2 || L[0] != "Speed:" || (speed = L[1].toDouble()) <= 0.0) error = true;
    Buffer = S.readLine();
    L = Buffer.split('\t');
    if (L.size() != 2 || L[0] != "StepSize:" || (stepSize = L[1].toDouble()) <= 0.0) error = true;
    Buffer = S.readLine();
    if (Buffer == "is moving") isMoving = true;
    else if (Buffer == "is not moving") isMoving = false;
    else error = true;
    if ((Buffer = S.readLine()) != " xp \t yp \t zp \t X \t Y \t Z \t vX \t vY \t vZ \t aaX \t aaY \t aaZ \t lX \t lY \t lZ \t lvX \t lvY \t lvZ ") error = true;
    Particle* part = Calc->getParticles(N);
    for (n=0; n<N && (!error); ++n)
    {
        L = S.readLine().split('\t');
        if (L.size() == 18)
        {
            part[n].R.setX(L[3].toDouble());
            part[n].R.setY(L[4].toDouble());
            part[n].R.setZ(L[5].toDouble());
            part[n].v.setX(L[6].toDouble());
            part[n].v.setY(L[7].toDouble());
            part[n].v.setZ(L[8].toDouble());
            part[n].aa.setX(L[9].toDouble());
            part[n].aa.setY(L[10].toDouble());
            part[n].aa.setZ(L[11].toDouble());
            part[n].lR.setX(L[12].toDouble());
            part[n].lR.setY(L[13].toDouble());
            part[n].lR.setZ(L[14].toDouble());
            part[n].lv.setX(L[15].toDouble());
            part[n].lv.setY(L[16].toDouble());
            part[n].lv.setZ(L[17].toDouble());
        }
        else error = true;
    }
    if (error)
    {
        QMessageBox::information(this, "Error reading file", "The format of the snapshot file is not supported.");
        if (n>0) Calc->initialize();
    }
    else
    {
        Calc->setStepSize(stepSize);
        Calc->setSpeed(speed);
        if (isMoving != Calc->getMove()) Calc->move();
        Calc->start();
    }
}

void Window::setLayerDistance(const double newDistance)
{
    if (nullptr == mNetworkClient) Calc->setLayerDistance(newDistance);
}

void Window::SetEnergyDefinitionAxis(const int particeIndex, const Vector &direction, Vector &end1, Vector &end2)
{
    Calc->CalcEndpointsOfEnergyDefinitionAxis(particeIndex, direction, end1, end2);
    mMarkedParticle = Calc->TranslateParticleIndex(mMarkedIndex = particeIndex);
    mAxisStart = end1;
    mAxisEnd = end2;
}

void Window::GetAxisEnergies(PotentialDefinerInputData &data)
{
    Calc->GetAxisEnergies(data);
}
