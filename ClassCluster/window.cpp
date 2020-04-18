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
#include <QMessageBox>


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
        S << (Calc->getMove() ? "is moving\n" : "is not moving\n");
        S << " xp \t yp \t zp \t X \t Y \t Z \t vX \t vY \t vZ \t aaX \t aaY \t aaZ \t lX \t lY \t lZ \t lvX \t lvY \t lvZ \n";
        for (int n=0; n<N; ++n)
            S << P[n].xp << "\t" << P[n].yp << "\t" << P[n].zp << "\t" << QString::number(P[n].X, 'f', 12) << "\t" << QString::number(P[n].Y, 'f', 12) << "\t" << QString::number(P[n].Z, 'f', 12)
              << "\t" << QString::number(P[n].vX, 'f', 12) << "\t" << QString::number(P[n].vY, 'f', 12) << "\t" << QString::number(P[n].vZ, 'f', 12)
              << "\t" << QString::number(P[n].aaX, 'f', 12) << "\t" << QString::number(P[n].aaY, 'f', 12) << "\t" << QString::number(P[n].aaZ, 'f', 12)
              << "\t" << QString::number(P[n].lX, 'f', 12) << "\t" << QString::number(P[n].lY, 'f', 12) << "\t" << QString::number(P[n].lZ, 'f', 12)
              << "\t" << QString::number(P[n].lvX, 'f', 12) << "\t" << QString::number(P[n].lvY, 'f', 12) << "\t" << QString::number(P[n].lvZ, 'f', 12) << "\n";
    }
}

void Window::restoreSnapShot(bool &isMoving)
{
    int n, N;
    double speed, stepSize;
    bool error = false;
    Calc->initialize();
    QString FileN = QFileDialog::getOpenFileName(this, "Select the snapshot to restore");
    if (FileN.isEmpty()) return;
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
            part[n].xp = L[0].toInt();
            part[n].yp = L[1].toInt();
            part[n].zp = L[2].toInt();
            part[n].X = L[3].toDouble();
            part[n].Y = L[4].toDouble();
            part[n].Z = L[5].toDouble();
            part[n].vX = L[6].toDouble();
            part[n].vY = L[7].toDouble();
            part[n].vZ = L[8].toDouble();
            part[n].aaX = L[9].toDouble();
            part[n].aaY = L[10].toDouble();
            part[n].aaZ = L[11].toDouble();
            part[n].lX = L[12].toDouble();
            part[n].lY = L[13].toDouble();
            part[n].lZ = L[14].toDouble();
            part[n].lvX = L[15].toDouble();
            part[n].lvY = L[16].toDouble();
            part[n].lvZ = L[17].toDouble();
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
