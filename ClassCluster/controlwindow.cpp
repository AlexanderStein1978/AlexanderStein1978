#include "controlwindow.h"

#include "window.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCloseEvent>


ControlWindow::ControlWindow() : window(new Window)
{
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(Start = new QPushButton("Start", this), 0, 0);
    L->addWidget(Restart = new QPushButton("Restart", this), 0, 1);
    L->addWidget(Rotate = new QPushButton("Rotate", this), 0, 2);
    L->addWidget(Move = new QPushButton("Move", this), 0, 3);
    L->addWidget(WriteSnapShot = new QPushButton("Write snapshot", this), 1, 0, 1, 2);
    L->addWidget(RestoreSnapShot = new QPushButton("Restore snapshot", this), 1, 2, 1, 2);
    L->addWidget(new QLabel("Speed:", this), 2, 0);
    L->addWidget(Speed = new QLineEdit(QString::number(1e3, 'f', 3), this), 2, 1);
    L->addWidget(new QLabel("Step size:", this), 2, 2);
    L->addWidget(StepE = new QLineEdit(QString::number(0.001, 'f', 3), this), 2, 3);
    L->addWidget(new QLabel("Energy:", this), 3, 0);
    L->addWidget(EnE = new QLineEdit(QString::number(window->getEnergy(), 'f', 3), this), 3, 1, 1, 3);
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
    if (NULL != window) delete window;
}

void ControlWindow::closeEvent(QCloseEvent *event)
{
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
        EnE->setText(QString::number(window->setEnergy(EnE->text().toDouble()), 'f', 3));
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
