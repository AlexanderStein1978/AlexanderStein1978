#include "window.h"
#include "Calculation.h"
#include "Picture.h"
#include "particle.h"
#include "vector.h"

#include <QCloseEvent>
#include <QPainter>
#include <QRect>
#include <QFileDialog>
#include <QTextStream>
#include <QBoxLayout>
#include <QMessageBox>


Window::Window(PotStruct *PotSs) : mPos(nullptr), mN(0)
{
    Calc = new Calculation(PotSs);
    QBoxLayout *L = new QBoxLayout(QBoxLayout::LeftToRight, this);
    L->addWidget(Pict = new Picture(this));
    int XSize, YSize;
    Calc->getSize(XSize, YSize);
    setMinimumSize(XSize, YSize);
    setMaximumSize(XSize, YSize);

    connect(Calc, SIGNAL(WriteSnapShot(Particle*, int)), this, SLOT(writeSnapShot(Particle*, int)));
    connect(Calc, SIGNAL(PictureChanged(Vector*, int)),
            this, SLOT(draw(Vector*, int)));
    connect(Calc, SIGNAL(EnergiesChanged(double,double)), this, SIGNAL(EnergiesChanged(double,double)));
}

Window::~Window()
{
    stopCalc();
    delete Calc;
    destroyData();
}

void Window::destroyData()
{
    if (0 != mN) delete[] mPos;
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

void Window::draw(Vector *Pos, int N)
{
    Calc->mutex.lock();
    if (mN != N)
    {
        destroyData();
        mN = N;
        mPos = new Vector[N];
    }
    memcpy(mPos, Pos, N * sizeof(Vector));
    Calc->mutex.unlock();
    update();
}

void Window::paintEvent(QPaintEvent *e)
{
    if (0 != mN)
    {
        int w = Pict->width(), he = Pict->height(), col, n, z, MaxZ;
        double ScF;
        Calc->getScales(ScF, MaxZ);
        double zSc = 255.0 / MaxZ, xc, yc;
        QPainter Paint(Pict->getPixmap());
        Paint.eraseRect(0, 0, w, he);
        for (n=0, z=-100; n < mN; n++)
        {
            if (int(mPos[n].Z() * zSc) != z)
            {
                col = 254 - (z = int(mPos[n].Z() * zSc));
                if (col < 0) col = 0;
                if (col > 254) col = 254;
                Paint.setPen(QColor(col, col, col));
                Paint.setBrush(QColor(col, col, col));
            }
            xc = mPos[n].X() * ScF;
            yc = mPos[n].Y() * ScF;
            Paint.drawEllipse(QRectF(xc - 5.0, yc - 5.0, 10.0, 10.0));
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

double Window::getPotentialEnergy() const
{
    return Calc->getPotentialEnergy();
}

double Window::getKineticEnergy() const
{
    return Calc->getKineticEnergy();
}

double Window::setKineticEnergy(const double T)
{
    return Calc->setKineticEnergy(T);
}

void Window::setParticleWatchPoint(WatchPoint *point)
{
    Calc->setParticleWatchPoint(point);
}

void Window::setParticleWatch(const int indexToWatch)
{
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
    Calc->setPotRangeScale(newScale);
}

void Window::setStepSize(const double size)
{
    Calc->setStepSize(size);
}

void Window::triggerSnapShot()
{
    Calc->triggerSnapShot();
}

void Window::setPotential(const Calculation::PotRole role, PotStruct &Pot)
{
    Calc->setPotential(role, Pot);
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
            S << P[n].xp << "\t" << P[n].yp << "\t" << P[n].zp << "\t" << QString::number(P[n].R.X(), 'f', 12) << "\t" << QString::number(P[n].R.Y(), 'f', 12) << "\t"
              << QString::number(P[n].R.Z(), 'f', 12) << "\t" << QString::number(P[n].v.X(), 'f', 12) << "\t" << QString::number(P[n].v.Y(), 'f', 12) << "\t"
              << QString::number(P[n].v.Z(), 'f', 12) << "\t" << QString::number(P[n].aa.X(), 'f', 12) << "\t" << QString::number(P[n].aa.Y(), 'f', 12) << "\t"
              << QString::number(P[n].aa.Z(), 'f', 12) << "\t" << QString::number(P[n].lR.X(), 'f', 12) << "\t" << QString::number(P[n].lR.Y(), 'f', 12) << "\t"
              << QString::number(P[n].lR.Z(), 'f', 12) << "\t" << QString::number(P[n].lv.X(), 'f', 12) << "\t" << QString::number(P[n].lv.Y(), 'f', 12) << "\t"
              << QString::number(P[n].lv.Z(), 'f', 12) << "\n";
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
    Calc->setLayerDistance(newDistance);
}
