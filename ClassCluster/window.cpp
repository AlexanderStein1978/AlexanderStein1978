#include "window.h"
#include "Calculation.h"
#include "picture.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCloseEvent>
#include <QPainter>
#include <QRect>


Window::Window()
{
    Calc = new Calculation;
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(Start = new QPushButton("Start", this), 0, 0);
    L->addWidget(Restart = new QPushButton("Restart", this), 0, 1);
    L->addWidget(Rotate = new QPushButton("Rotate", this), 0, 2);
    L->addWidget(Move = new QPushButton("Move", this), 0, 3);
    L->addWidget(new QLabel("Speed:", this), 0, 4);
    L->addWidget(Speed = new QLineEdit(QString::number(Calc->getSpeed(), 'f', 3), this), 0, 5);
    L->addWidget(new QLabel("Step size:", this), 0, 6);
    L->addWidget(StepE = new QLineEdit(QString::number(Calc->getStepSize(), 'f', 3), this), 0, 7);
    L->addWidget(new QLabel("Energy:", this), 0, 8);
    L->addWidget(EnE = new QLineEdit(QString::number(Calc->getEnergy(), 'f', 3), this), 0, 9);
    L->addWidget(End = new QPushButton("Quit", this), 0, 10);
    L->addWidget(Pict = new Picture(this), 1, 0, 1, 11);
    L->setRowStretch(1, 1);
    int XSize, YSize;
    Calc->getSize(XSize, YSize);
    setMinimumSize(XSize, YSize);
    setMaximumSize(XSize, YSize);
    connect(Start, SIGNAL(clicked()), this, SLOT(run()));
    connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));
    connect(Rotate, SIGNAL(clicked()), Calc, SLOT(rotate()));
    connect(Move, SIGNAL(clicked()), this, SLOT(move()));
    connect(Speed, SIGNAL(editingFinished()), this, SLOT(speedChanged()));
    connect(End, SIGNAL(clicked()), this, SLOT(close()));
    connect(Calc, SIGNAL(PictureChanged(double*, double*, double*, int)),
            this, SLOT(draw(double*, double*, double*, int)));
}

Window::~Window()
{
    stopCalc();
    delete Calc;
}

void Window::stopCalc()
{
    if (Calc->isRunning())
    {
        Calc->stop();
        Calc->wait();
    }
}

void Window::closeEvent(QCloseEvent* event)
{
    stopCalc();
    event->accept();
}

void Window::draw(double* XP, double* YP, double* ZP, int N)
{
    int w = Pict->width(), he = Pict->height(), col, n, z, MaxZ;
    double ScF;
    Calc->getScales(ScF, MaxZ);
    double zSc = 255.0 / MaxZ, xc, yc;
    QPainter Paint(Pict->getPixmap());
    Paint.eraseRect(0, 0, w, he);
    Calc->mutex.lock();
    for (n=0, z=-100; n<N; n++)
    {
        if (int(ZP[n] * zSc) != z)
        {
            col = 254 - (z = int(ZP[n] * zSc));
            if (col < 0) col = 0;
            if (col > 254) col = 254;
            Paint.setPen(QColor(col, col, col));
            Paint.setBrush(QColor(col, col, col));
        }
        xc = XP[n] * ScF;
        yc = YP[n] * ScF;
        Paint.drawEllipse(QRectF(xc - 5.0, yc - 5.0, 10.0, 10.0));
    }
    Calc->mutex.unlock();
    Pict->repaint();
}

void Window::run()
{
    if (Calc->isRunning())
    {
        Calc->stop();
        Start->setText("Start");
    }
    else
    {
        EnE->setText(QString::number(Calc->setEnergy(EnE->text().toDouble()), 'f', 3));
        Calc->setStepSize(StepE->text().toDouble());
        Calc->start();
        Start->setText("Stop");
    }
}

void Window::restart()
{
    stopCalc();
    Calc->initialize();
    run();
}

void Window::move()
{
    if (Calc->getMove()) Move->setText("Move");
    else Move->setText("Hold");
    Calc->move();
}

void Window::speedChanged()
{
    Calc->setSpeed(Speed->text().toDouble());
}
