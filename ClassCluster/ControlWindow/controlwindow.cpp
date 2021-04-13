#include "controlwindow.h"
#include "Calculation.h"
#include "window.h"
#include "potcontrol.h"
#include "potstruct.h"
#include "potentialplot.h"
#include "MainWindow.h"
#include "particlewatchtable.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCloseEvent>
#include <QTextStream>
#include <QDir>


ControlWindow::ControlWindow(MainWindow * const mw) : window(nullptr), TEdit(new QLineEdit(this)),
    PotControls(new PotControl*[Calculation::NumPot]), Plot(nullptr), SettingsFileName("../../../Physics/ClassCluster/Data/Settings.dat"),
    MW(mw), ProgramPath(QDir::currentPath())
{
    QFile settingsFile(SettingsFileName);
    QString speed(QString::number(1e3, 'f', 3)), stepSize(QString::number(1e-3, 'f', 3)), kineticEnergy("-1.0"), rangeScale(QString::number(1.0, 'f', 3)), layerDistance("5.657");
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n] = new PotControl(this, mw);
    if (settingsFile.exists())
    {
        settingsFile.open(QIODevice::ReadOnly);
        QTextStream S(&settingsFile);
        QStringList settings = S.readLine().split('\t');
        if (settings.size() >= 4u)
        {
            speed = settings[0];
            stepSize = settings[1];
            kineticEnergy = settings[2];
            if (kineticEnergy.toDouble() < 0.0) kineticEnergy = "0.0";
            TEdit->setText(kineticEnergy);
            rangeScale = settings[3];
            if (settings.size() >= 5u) layerDistance = settings[4];
            for (int n=0; n < Calculation::NumPot && !S.atEnd(); ++n)
            {
                const QString data = S.readLine();
                PotControls[n]->Init(data);
            }
        }
    }
    prepareWindow();
    const double cEnergy = window->getPotentialEnergy() + kineticEnergy.toDouble();
    QGridLayout *L = new QGridLayout(this), *SettingsLayout = new QGridLayout, *PotLayout = new QGridLayout;
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
    SettingsLayout->addWidget(Speed = new QLineEdit(speed, this), 0, 1);
    SettingsLayout->addWidget(new QLabel("Step size:", this), 0, 2);
    SettingsLayout->addWidget(StepE = new QLineEdit(stepSize, this), 0, 3);
    SettingsLayout->addWidget(new QLabel("Potential range:"), 0, 4);
    SettingsLayout->addWidget(PotRangeScaleEdit = new QLineEdit(rangeScale, this), 0, 5);
    SettingsLayout->addWidget(new QLabel("Layer distance:", this), 1, 0);
    SettingsLayout->addWidget(LayerDistanceEdit = new QLineEdit(layerDistance, this), 1, 1);
    SettingsLayout->addWidget(new QLabel("Kintetic energy:", this), 1, 2);
    SettingsLayout->addWidget(TEdit, 1, 3);
    SettingsLayout->addWidget(new QLabel("Total energy:", this), 1, 4);
    SettingsLayout->addWidget(EnE = new QLineEdit(QString::number(cEnergy, 'g', 6), this), 1, 5);
    SettingsLayout->addWidget(PotentialEnergyLabel = new QLabel("Potential energy: ", this), 2, 0, 1, 2);
    SettingsLayout->addWidget(KineticEnergyLabel = new QLabel("Kinetic energy: ", this), 2, 2, 1, 2);
    SettingsLayout->addWidget(TotalEnergyLabel = new QLabel("Total energy: ", this), 2, 4, 1, 2);
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
    for (int n=0; n < Calculation::NumPot; ++n) connect(PotControls[n], SIGNAL(Change()), this, SLOT(EnergyRelevantValueChanged()));
    setFocusPolicy(Qt::StrongFocus);
}

ControlWindow::~ControlWindow()
{
    for (int n=0; n < Calculation::NumPot; ++n) delete PotControls[n];
    delete[] PotControls;
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
    if (window->isRunning())
    {
        window->stop();
        Start->setText("Start");
    }
    else
    {
        if (!window->isVisible()) window->show();
        start();
        Start->setText("Stop");
    }
}

void ControlWindow::start()
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
    window->start();
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
    S << Speed->text() << '\t' << StepE->text() << '\t' << TEdit->text() << '\t' << PotRangeScaleEdit->text() << '\t' << LayerDistanceEdit->text() << '\n';
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->Serialize(S, ProgramPath);
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
