#include "controlwindow.h"
#include "Calculation.h"
#include "window.h"
#include "potcontrol.h"
#include "potstruct.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCloseEvent>
#include <QTextStream>


ControlWindow::ControlWindow() : window(nullptr), PotControls(new PotControl*[Calculation::NumPot]), SettingsFileName("../../../Physics/ClassCluster/Settings.dat")
{
    QFile settingsFile(SettingsFileName);
    QString speed(QString::number(1e3, 'f', 3)), stepSize(QString::number(1e-3, 'f', 3)), energy(QString::number(window->getEnergy(), 'g', 3)), rangeScale(QString::number(1.0, 'f', 3));
    PotStruct PotSs[Calculation::NumPot];
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n] = new PotControl(this);
    if (settingsFile.exists())
    {
        settingsFile.open(QIODevice::ReadOnly);
        QTextStream S(&settingsFile);
        QStringList settings = S.readLine().split('\t', QString::SkipEmptyParts);
        if (settings.size() == 4u)
        {
            speed = settings[0];
            stepSize = settings[1];
            energy = settings[2];
            rangeScale = settings[3];
            for (int n=0; n < Calculation::NumPot && !S.atEnd(); ++n)
            {
                const QString data = S.readLine();
                PotControls[n]->Init(data);
                PotControls[n]->FillStruct(PotSs[n]);
            }
        }
    }
    window = new Window(PotSs);
    QGridLayout *L = new QGridLayout(this), *PotLayout = new QGridLayout;
    L->setColumnStretch(0, 1);
    L->setColumnStretch(1, 1);
    L->setColumnStretch(2, 1);
    L->setColumnStretch(3, 1);
    L->addWidget(Start = new QPushButton("Start", this), 0, 0);
    L->addWidget(Restart = new QPushButton("Restart", this), 0, 1);
    L->addWidget(Rotate = new QPushButton("Rotate", this), 0, 2);
    L->addWidget(Move = new QPushButton("Move", this), 0, 3);
    L->addWidget(WriteSnapShot = new QPushButton("Write snapshot", this), 1, 0, 1, 2);
    L->addWidget(RestoreSnapShot = new QPushButton("Restore snapshot", this), 1, 2, 1, 2);
    L->addWidget(new QLabel("Speed:", this), 2, 0);
    L->addWidget(Speed = new QLineEdit(speed, this), 2, 1);
    L->addWidget(new QLabel("Step size:", this), 2, 2);
    L->addWidget(StepE = new QLineEdit(stepSize, this), 2, 3);
    L->addWidget(new QLabel("Energy:", this), 3, 0);
    L->addWidget(EnE = new QLineEdit(energy, this), 3, 1);
    L->addWidget(new QLabel("Potential range:"), 3, 2);
    L->addWidget(PotRangeScaleEdit = new QLineEdit(rangeScale, this), 3, 3);
    L->setRowMinimumHeight(4, 20);
    L->addWidget(new QLabel("Interaction potentials:", this), 5, 0, 1, 4);
    L->addLayout(PotLayout, 6, 0, 1, 4);
    PotLayout->addWidget(new QLabel("Closest 2", this), 0, 0);
    PotLayout->addWidget(new QLabel("Next 2", this), 0, 0);
    PotLayout->addWidget(new QLabel("Remaining", this), 0, 0);
    for (int n=0; n < Calculation::NumPot; ++n) PotControls[n]->FillLayout(PotLayout, n);
    connect(Start, SIGNAL(clicked()), this, SLOT(run()));
    connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));
    connect(Rotate, SIGNAL(clicked()), this, SLOT(rotate()));
    connect(Move, SIGNAL(clicked()), this, SLOT(move()));
    connect(Speed, SIGNAL(editingFinished()), this, SLOT(speedChanged()));
    connect(WriteSnapShot, SIGNAL(clicked()), this, SLOT(writeSnapShot()));
    connect(RestoreSnapShot, SIGNAL(clicked()), this, SLOT(restoreSnapShot()));
}

ControlWindow::~ControlWindow()
{
    QFile file(SettingsFileName);
    file.open(QIODevice::WriteOnly);
    QTextStream S(&file);
    S << Speed->text() << '\t' << StepE->text() << '\t' << EnE->text() << '\t' << PotRangeScaleEdit << '\n';
    if (NULL != window) delete window;
    for (int n=0; n < Calculation::NumPot; ++n)
    {
        PotControls[n]->Serialize(S);
        delete PotControls[n];
    }
    delete[] PotControls;
}

void ControlWindow::closeEvent(QCloseEvent *event)
{
    for (int n=0; n < Calculation::NumPot; ++n) if (!PotControls[n]->canPotBeClosed())
    {
        event->ignore();
        return;
    }
    if (NULL != window) window->close();
    event->accept();
}

void ControlWindow::run()
{
    if (window->isRunning())
    {
        window->stop();
        Start->setText("Start");
    }
    else
    {
        if (!window->isVisible()) window->show();
        window->setPotentialRangeScale(PotRangeScaleEdit->text().toDouble());
        for (int n=0; n < Calculation::NumPot; ++n) if (PotControls[n]->isChanged())
        {
            PotStruct Pots;
            PotControls[n]->FillStruct(Pots);
            window->setPotential(static_cast<Calculation::PotRole>(n), Pots);
        }
        EnE->setText(QString::number(window->setEnergy(EnE->text().toDouble()), 'g', 3));
        window->setStepSize(StepE->text().toDouble());
        window->start();
        Start->setText("Stop");
    }
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
