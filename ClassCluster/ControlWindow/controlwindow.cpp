#include "controlwindow.h"
#include "Calculation.h"
#include "window.h"
#include "potcontrol.h"
#include "potstruct.h"
#include "potentialplot.h"
#include "MainWindow.h"
#include "potential.h"
#include "particlewatchtable.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCloseEvent>
#include <QTextStream>
#include <QDir>


ControlWindow::ControlWindow(MainWindow * const mw) : window(nullptr), StepE(new QLineEdit(QString::number(1e-3, 'f', 3), this)), EnE(new QLineEdit(this)), TEdit(new QLineEdit("-1.0", this)),
    Speed(new QLineEdit(QString::number(1e3, 'f', 3), this)), PotRangeScaleEdit(new QLineEdit(QString::number(1.0, 'f', 3), this)), LayerDistanceEdit(new QLineEdit("5.657", this)),
    IpAddressEdit(new QLineEdit("192.168.1.1", this)), ConnectionStatus(new QLabel("disconnected", this)), NetworkSelection(new QComboBox(this)), Connect(new QPushButton("Connect", this)),
    GetSettings(new QPushButton("Get settings", this)), PotControls(new PotControl*[Calculation::NumPot]), Plot(nullptr), MW(mw), SettingsFileName("../../../Physics/ClassCluster/Data/Settings.dat"), ProgramPath(QDir::currentPath())
{
    QFile settingsFile(SettingsFileName);
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n] = new PotControl(this, mw);
    if (settingsFile.exists())
    {
        settingsFile.open(QIODevice::ReadOnly);
        QTextStream S(&settingsFile);
        Init(S);
    }
    prepareWindow();
    const double cEnergy = window->getPotentialEnergy() + TEdit->text().toDouble();
    QGridLayout *L = new QGridLayout(this), *SettingsLayout = new QGridLayout, *PotLayout = new QGridLayout, *CornerLayout = new QGridLayout;
    L->setColumnStretch(0, 1);
    L->setColumnStretch(1, 1);
    L->setColumnStretch(2, 1);
    L->setColumnStretch(3, 1);
    L->addWidget(Start = new QPushButton("Start", this), 0, 0);
    L->addWidget(Restart = new QPushButton("Restart", this), 0, 1);
    L->addWidget(Rotate = new QPushButton("Rotate", this), 0, 2);
    L->addWidget(Move = new QPushButton("Move", this), 0, 3);
    L->addWidget(WriteSnapShot = new QPushButton("Write snapshot", this), 1, 0);
    L->addWidget(RestoreSnapShot = new QPushButton("Restore snapshot", this), 1, 1);
    L->addWidget(ShowParticleWatchWindow = new QPushButton("Show particle amplification watch window", this), 1, 2, 1, 2);
    L->addLayout(SettingsLayout, 2, 0, 1, 4);
    SettingsLayout->addWidget(new QLabel("Speed:", this), 0, 0);
    SettingsLayout->addWidget(Speed , 0, 1);
    SettingsLayout->addWidget(new QLabel("Step size:", this), 0, 2);
    SettingsLayout->addWidget(StepE, 0, 3);
    SettingsLayout->addWidget(new QLabel("Potential range:"), 0, 4);
    SettingsLayout->addWidget(PotRangeScaleEdit, 0, 5);
    SettingsLayout->addWidget(new QLabel("Layer distance:", this), 1, 0);
    SettingsLayout->addWidget(LayerDistanceEdit, 1, 1);
    SettingsLayout->addWidget(new QLabel("Kintetic energy:", this), 1, 2);
    SettingsLayout->addWidget(TEdit, 1, 3);
    SettingsLayout->addWidget(new QLabel("Total energy:", this), 1, 4);
    SettingsLayout->addWidget(EnE, 1, 5);
    EnE->setText(QString::number(cEnergy, 'g', 6));
    SettingsLayout->addWidget(PotentialEnergyLabel = new QLabel("Potential energy: ", this), 2, 0, 1, 2);
    SettingsLayout->addWidget(KineticEnergyLabel = new QLabel("Kinetic energy: ", this), 2, 2, 1, 2);
    SettingsLayout->addWidget(TotalEnergyLabel = new QLabel("Total energy: ", this), 2, 4, 1, 2);
    SettingsLayout->addWidget(NetworkSelection, 3, 0, 1, 2);
    NetworkSelection->addItems(QStringList() << "No network, use local calculation" << "Listen as calculation server" << "Connect to calculation server");
    SettingsLayout->addWidget(new QLabel("IP Address:", this), 3, 2);
    SettingsLayout->addWidget(IpAddressEdit, 3, 3);
    SettingsLayout->addWidget(Connect, 3, 4);
    Connect->setEnabled(false);
    CornerLayout->addWidget(ConnectionStatus, 0, 0);
    CornerLayout->addWidget(GetSettings, 0, 1);
    SettingsLayout->addLayout(CornerLayout, 3, 5);
    disconnected();
    L->setRowMinimumHeight(3, 20);
    L->addWidget(new QLabel("Interaction potentials:", this), 4, 0, 1, 4);
    L->addLayout(PotLayout, 5, 0, 1, 4);
    PotLayout->setColumnMinimumWidth(0, 60);
    PotLayout->setColumnStretch(1, 1);
    PotLayout->setColumnStretch(2, 10);
    PotLayout->setColumnStretch(3, 1);
    PotLayout->setColumnMinimumWidth(4, 20);
    PotLayout->setColumnStretch(5, 1);
    PotLayout->setColumnMinimumWidth(6, 20);
    PotLayout->setColumnStretch(7, 1);
    PotLayout->setColumnStretch(8, 1);
    PotLayout->addWidget(new QLabel("Closest 2:", this), Calculation::ClosestTwo, 0);
    PotLayout->addWidget(new QLabel("Next 2:", this), Calculation::NextTwo, 0);
    PotLayout->addWidget(new QLabel("Second order:", this), Calculation::SecondOrder, 0);
    PotLayout->addWidget(new QLabel("Remaining:", this), Calculation::Remaining, 0);
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->FillLayout(PotLayout, n);
    connect(Start, SIGNAL(clicked()), this, SLOT(run()));
    connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));
    connect(Rotate, SIGNAL(clicked()), this, SLOT(rotate()));
    connect(Move, SIGNAL(clicked()), this, SLOT(move()));
    connect(Speed, SIGNAL(editingFinished()), this, SLOT(speedChanged()));
    connect(WriteSnapShot, SIGNAL(clicked()), this, SLOT(writeSnapShot()));
    connect(RestoreSnapShot, SIGNAL(clicked()), this, SLOT(restoreSnapShot()));
    connect(ShowParticleWatchWindow, SIGNAL(clicked()), this, SLOT(showParticleWatchWindow()));
    connect(StepE, SIGNAL(editingFinished()), this, SLOT(ValueChanged()));
    connect(EnE, SIGNAL(editingFinished()), this, SLOT(EChanged()));
    connect(TEdit, SIGNAL(editingFinished()), this, SLOT(TChanged()));
    connect(PotRangeScaleEdit, SIGNAL(editingFinished()), this, SLOT(EnergyRelevantValueChanged()));
    connect(LayerDistanceEdit, SIGNAL(editingFinished()), this, SLOT(EnergyRelevantValueChanged()));
    connect(MW, SIGNAL(MainWindowCloses()), this, SLOT(saveSettings()));
    connect(window, SIGNAL(EnergiesChanged(double,double)), this, SLOT(UpdateEnergies(double,double)));
    connect(NetworkSelection, SIGNAL(currentIndexChanged(int)), this, SLOT(networkSelectionChanged(int)));
    connect(Connect, SIGNAL(clicked()), this, SLOT(connectToServer()));
    connect(GetSettings, SIGNAL(clicked()), this, SLOT(sendGetSettingsRequest()));
    connect(window, SIGNAL(IsConnectedToServer()), this, SLOT(connectionEstablished()));
    connect(window, SIGNAL(ConnectionFailed()), this, SLOT(disconnected()));
    connect(window, SIGNAL(IsRunning(bool)), this, SLOT(setIsRunning(bool)));
    connect(window, SIGNAL(SendSettings()), this, SLOT(getSettings()));
    connect(window, SIGNAL(ReceivedSetting(const QByteArray&)), this, SLOT(setSettings(const QByteArray&)));
    connect(window, SIGNAL(ReceivedPotential(const QByteArray&)), this, SLOT(setPotentialData(const QByteArray&)));
    connect(window, SIGNAL(ReceivedStartCommand()), this, SLOT(run()));
    for (int n=0; n < Calculation::NumPot; ++n) connect(PotControls[n], SIGNAL(Change()), this, SLOT(EnergyRelevantValueChanged()));
    setFocusPolicy(Qt::StrongFocus);
}

ControlWindow::~ControlWindow()
{
    for (int n=0; n < Calculation::NumPot; ++n) delete PotControls[n];
    delete[] PotControls;
}

void ControlWindow::Init(QTextStream& inStream)
{
    QStringList settings = inStream.readLine().split('\t');
    if (settings.size() >= 4)
    {
        Speed->setText(settings[0]);
        StepE->setText(settings[1]);
        QString kineticEnergy = settings[2];
        if (kineticEnergy.toDouble() < 0.0) kineticEnergy = "0.0";
        TEdit->setText(kineticEnergy);
        PotRangeScaleEdit->setText(settings[3]);
        if (settings.size() >= 5) LayerDistanceEdit->setText(settings[4]);
        for (int n=0; n < Calculation::NumPot && !inStream.atEnd(); ++n)
        {
            const QString data = inStream.readLine();
            PotControls[n]->Init(data);
        }
    }
}

void ControlWindow::Init(QString& data)
{
    QTextStream instream(&data, QIODevice::ReadOnly);
    Init(instream);
}

void ControlWindow::setSettings(const QByteArray& data)
{
    QTextStream stream(data, QIODevice::ReadOnly);
    Init(stream);
    window->setSpeed(Speed->text().toDouble());
    if (window->isRunning())
    {
        window->stop();
        start();
    }
}

void ControlWindow::setPotentialData(const QByteArray& data)
{
    QTextStream stream(data, QIODevice::ReadOnly);
    stream.readLine();
    QStringList NameL = stream.readLine().split(": ");
    if (NameL.size() >= 2)
    {
        QString Name = NameL[1];
        if (!Name.isEmpty())
        {
            Potential* pot = MW->getPotential(Name);
            if (nullptr == pot)
            {
                pot = MW->CreatePotential();
                if (nullptr == pot) return;
                pot->show();
            }
            stream.seek(0);
            pot->init(stream);
        }
    }
}

void ControlWindow::networkSelectionChanged(int index)
{
    switch(index)
    {
    case NetSelLocal:
        IpAddressEdit->setEnabled(false);
        Connect->setEnabled(false);
        window->switchBackToLocalCalulations();
        window->stopListeningAsCalculationServer();
        break;
    case NetSelServer:
        IpAddressEdit->setEnabled(true);
        Connect->setEnabled(true);
        Connect->setText("Listen");
        window->switchBackToLocalCalulations();
        break;
    case NetSelClient:
        IpAddressEdit->setEnabled(true);
        Connect->setEnabled(true);
        Connect->setText("Connect");
        window->stopListeningAsCalculationServer();
        break;
    }
    setConnectionStatus(false);
}

void ControlWindow::connectToServer()
{
    if (NetworkSelection->currentIndex() == 2) window->connectToCalculationServer(IpAddressEdit->text());
    else setConnectionStatus(window->listenAsCalculationServer(IpAddressEdit->text()));
}

void ControlWindow::setConnectionStatus(bool connected)
{
    QPalette palette(ConnectionStatus->palette());
    if (connected) palette.setColor(QPalette::WindowText, QColor(0, 255, 0));
    else palette.setColor(QPalette::WindowText, QColor(255, 0, 0));
    ConnectionStatus->setPalette(palette);
    if (NetworkSelection->currentIndex() == 2) ConnectionStatus->setText(connected ? "connected" : "disconnected");
    else ConnectionStatus->setText(connected ? "listening" : "failed");
}

void ControlWindow::UpdateEnergies(double kineticEnergy, double totalEnergy)
{
    PotentialEnergyLabel->setText("Potential energy: " + QString::number(totalEnergy - kineticEnergy, 'g', 3));
    KineticEnergyLabel->setText("Kinetic energy: " + QString::number(kineticEnergy, 'g', 3));
    TotalEnergyLabel->setText("Total energy: " + QString::number(totalEnergy, 'g', 3));
}

void ControlWindow::focusInEvent(QFocusEvent *event)
{
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->UpdatePotentialBox();
    event->accept();
}

double ControlWindow::getRe() const
{
    return window->getRe();
}

void ControlWindow::prepareWindow()
{
    PotStruct PotSs[Calculation::NumPot];
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->FillStruct(PotSs[n]);
    window = new Window(PotSs);
    MW->showMDIChild(window);
}

void ControlWindow::run()
{
    if (nullptr == window) prepareWindow();
    if (window->isRunning()) window->stop();
    else
    {
        if (!window->isVisible()) window->show();
        start();
    }
}

void ControlWindow::setIsRunning(bool isRunning)
{
    if (isRunning) Start->setText("Stop");
    else Start->setText("Start");
}

void ControlWindow::start()
{
    if (NetworkSelection->currentIndex() == NetSelClient) getSettings();
    else
    {
        window->setPotentialRangeScale(PotRangeScaleEdit->text().toDouble());
        window->setLayerDistance(LayerDistanceEdit->text().toDouble());
        for (int n=0; n < Calculation::NumPot; ++n) if (PotControls[n]->isChangedSinceLastRun())
        {
            PotStruct Pots;
            PotControls[n]->FillStruct(Pots);
            window->setPotential(static_cast<Calculation::PotRole>(n), Pots);
        }
        EnE->setText(QString::number(window->setKineticEnergy(TEdit->text().toDouble()), 'g', 3));
        window->setStepSize(StepE->text().toDouble());
    }
    window->start();
}

void ControlWindow::getSettings()
{
    QByteArray data;
    QTextStream stream(&data);
    Serialize(stream);
    stream.flush();
    window->SendSettings(data);
    for (int n=0; n < Calculation::NumPot; ++n) if (PotControls[n]->isChangedSinceLastRun())
    {
        data.clear();
        Potential* pot = PotControls[n]->getPotential();
        if (nullptr != pot)
        {
            pot->serialize(stream);
            stream.flush();
            window->SendPotential(data);
        }
    }
}

void ControlWindow::sendGetSettingsRequest()
{
    window->SendGetSettingsRequest();
}

void ControlWindow::restart()
{
    if (NULL != window) window->reset();
    run();
}

void ControlWindow::move()
{
    if (window->isMoving()) Move->setText("Move");
    else Move->setText("Hold");
    window->move();
}

void ControlWindow::rotate()
{
    if (NULL == window) run();
    window->rotate();
}

void ControlWindow::showPotential(Potential * const pot, const bool plot)
{
    if (nullptr == Plot)
    {
        Plot = new PotentialPlot;
        Plot->setShowHistory(false);
        Plot->setShowPoints(true);
        MW->showMDIChild(Plot);
        connect(Plot, SIGNAL(closing()), this, SLOT(plotClosing()));
    }
    else if (!Plot->isVisible()) Plot->show();
    if (plot) Plot->plotPotential(pot);
    else Plot->removePotential(pot);
}

void ControlWindow::plotClosing()
{
    disconnect(Plot, SIGNAL(closing()), this, SLOT(plotClosing()));
    Plot = nullptr;
    for (int i=0; i < Calculation::NumPot; ++i) PotControls[i]->PLotCloses();
}

void ControlWindow::showParticleWatchWindow()
{
    MW->showMDIChild(new ParticleWatchTable(window, MW));
}

void ControlWindow::speedChanged()
{
    if (NULL == window) run();
    window->setSpeed(Speed->text().toDouble());
}

void ControlWindow::writeSnapShot()
{
    if (NULL == window) run();
    window->triggerSnapShot();
}

void ControlWindow::restoreSnapShot()
{
    bool isMoving;
    if (NULL != window)
    {
        if (window->isRunning()) window->stopCalc();
    }
    else window = new Window;
    if (!window->isVisible()) window->show();
    window->restoreSnapShot(isMoving);
    Start->setText("Stop");
    if (isMoving) Move->setText("Hold");
    else Move->setText("Move");
}

void ControlWindow::saveSettings()
{
    QFile file(SettingsFileName);
    file.open(QIODevice::WriteOnly);
    QTextStream S(&file);
    Serialize(S);
}

void ControlWindow::Serialize(QTextStream& outStream)
{
    outStream << Speed->text() << '\t' << StepE->text() << '\t' << TEdit->text() << '\t' << PotRangeScaleEdit->text() << '\t' << LayerDistanceEdit->text() << '\n';
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->Serialize(outStream, ProgramPath);
}

void ControlWindow::EChanged()
{
    bool wasWindowRunning(stopIfItsRunning());
    double V = window->getPotentialEnergy(), T(EnE->text().toDouble() - V);
    if (T < 0.0)
    {
        T = 0.0;
        EnE->setText(QString::number(T));
    }
    TEdit->setText(QString::number(T));
    if (wasWindowRunning) start();
}

void ControlWindow::TChanged()
{
    EnergyRelevantValueChanged();
}

void ControlWindow::ValueChanged()
{
    if (window->isRunning())
    {
        window->stopCalc();
        start();
    }
}

void ControlWindow::EnergyRelevantValueChanged()
{
    bool wasWindowRunning(stopIfItsRunning());
    window->setLayerDistance(LayerDistanceEdit->text().toDouble());
    EnE->setText(QString::number(TEdit->text().toDouble() + window->getPotentialEnergy()));
    if (wasWindowRunning) start();
}

bool ControlWindow::stopIfItsRunning()
{
    if (window->isRunning())
    {
        window->stopCalc();
        return true;
    }
    return false;
}
