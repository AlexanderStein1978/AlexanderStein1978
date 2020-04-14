#include "window.h"
#include "Calculation.h"
#include "picture.h"
#include "particle.h"

#include <QCloseEvent>
#include <QPainter>
#include <QRect>
#include <QFileDialog>
#include <QTextStream>
#include <QBoxLayout>


Window::Window()
{
    Calc = new Calculation;
    QBoxLayout *L = new QBoxLayout(QBoxLayout::LeftToRight, this);
    L->addWidget(Pict = new Picture(this));
    int XSize, YSize;
    Calc->getSize(XSize, YSize);
    setMinimumSize(XSize, YSize);
    setMaximumSize(XSize, YSize);

    connect(Calc, SIGNAL(WriteSnapShot(Particle*, int)), this, SLOT(writeSnapShot(Particle*, int)));
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

void Window::start()
{
    Calc->start();
}

void Window::stop()
{
    Calc->stop();
}

bool Window::isRunning() const
{
    return Calc->isRunning();
}

void Window::reset()
{
    stopCalc();
    Calc->initialize();
}

void Window::rotate()
{
    Calc->rotate();
}

void Window::move()
{
    Calc->move();
}

bool Window::isMoving() const
{
    return Calc->getMove();
}

void Window::setSpeed(const double newSpeed)
{
    Calc->setSpeed(newSpeed);
}

double Window::setEnergy(const double E)
{
    return Calc->setEnergy(E);
}

void Window::setStepSize(const double size)
{
    Calc->setStepSize(size);
}

void Window::triggerSnapShot()
{
    Calc->triggerSnapShot();
}

void Window::writeSnapShot(Particle *P, int N)
{
    QString FN = QFileDialog::getSaveFileName(this, "Select filename and path for snapshot");
    if (!FN.isEmpty())
    {
        QFile file(FN);
        file.open(QIODevice::WriteOnly);
        QTextStream S(&file);
        S << "Speed:\t" << QString::number(Calc->getSpeed(), 'f', 12) << "\n";
        S << "StepSize:\t" << QString::number(Calc->getStepSize(), 'f', 12) << "\n";
        S << " xp \t yp \t zp \t X \t Y \t Z \t vX \t vY \t vZ \t aaX \t aaY \t aaZ \t lX \t lY \t lZ \t lvX \t lvY \t lvZ \n";
        for (int n=0; n<N; ++n)
            S << P[n].xp << "\t" << P[n].yp << "\t" << P[n].zp << "\t" << QString::number(P[n].X, 'f', 12) << "\t" << QString::number(P[n].Y, 'f', 12) << "\t" << QString::number(P[n].Z, 'f', 12)
              << "\t" << QString::number(P[n].vX, 'f', 12) << "\t" << QString::number(P[n].vY, 'f', 12) << "\t" << QString::number(P[n].vZ, 'f', 12)
              << "\t" << QString::number(P[n].aaX, 'f', 12) << "\t" << QString::number(P[n].aaY, 'f', 12) << "\t" << QString::number(P[n].aaZ, 'f', 12)
              << "\t" << QString::number(P[n].lX, 'f', 12) << "\t" << QString::number(P[n].lY, 'f', 12) << "\t" << QString::number(P[n].lZ, 'f', 12)
              << "\t" << QString::number(P[n].lvX, 'f', 12) << "\t" << QString::number(P[n].lvY, 'f', 12) << "\t" << QString::number(P[n].lvZ, 'f', 12) << "\n";
    }
}
