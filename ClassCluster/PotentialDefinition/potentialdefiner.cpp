//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "potentialdefiner.h"
#include "window.h"

#include <QLineEdit>
#include <QPainter>
#include <QGridLayout>
#include <QLabel>


PotentialDefiner::PotentialDefiner(Window* window, MainWindow* MW) : DiagWindow(SimpleDiagWindow, MW, "", "*.*", 1), mWindow(window), mCalcB(new QPushButton("Calculate", this)),
    mParticleIndexInput(new QLineEdit("0", this)), mDirectionXInput(new QLineEdit("0.0", this)), mDirectionYInput(new QLineEdit("0.0", this)), mDirectionZInput(new QLineEdit("0.0", this))
{
    setWindowTitle("Potential definition window");
    xStartLabel->setText("R [A]:    from ");
    xStart->setText("0");
    xStopLabel->setText(" to ");
    xStop->setText("100");
    yStartLabel->setText("Energie [cm-1]  :    from ");
    yStart->setText("0");
    yStopLabel->setText(" to ");
    yStop->setText("5000");
    XUnit = "% of axis length";
    YUnit = "Energy [cm^{-1}]";
    XMin = 0.0;
    XMax = 100.0;
    YMin = -10000.0;
    YMax = 50000.0;
    mParticleIndexInput->setValidator(new QIntValidator(0, window->getNumParticles() - 1, mParticleIndexInput));
    mDirectionXInput->setValidator(new QDoubleValidator(mDirectionXInput));
    mDirectionYInput->setValidator(new QDoubleValidator(mDirectionYInput));
    mDirectionZInput->setValidator(new QDoubleValidator(mDirectionZInput));
    QGridLayout *layout = new QGridLayout;
    layout->addWidget(new QLabel("Particle:", this), 0, 0);
    layout->addWidget(mParticleIndexInput, 0, 1);
    layout->addWidget(new QLabel("Direction x:"), 0, 2);
    layout->addWidget(mDirectionXInput, 0, 3);
    layout->addWidget(new QLabel("y:", this), 0, 4);
    layout->addWidget(mDirectionYInput, 0, 5);
    layout->addWidget(new QLabel("z:", this), 0, 6);
    layout->addWidget(mDirectionZInput, 0, 7);
    layout->addWidget(mCalcB, 0, 8);
    SpektrumLayout->addLayout(layout, 0, 0, 1, 9);
    connect(mCalcB, SIGNAL(clicked()), this, SLOT(Calculate()));
}

void PotentialDefiner::Calculate()
{
    Vector direction(mDirectionXInput->text().toDouble(), mDirectionYInput->text().toDouble(), mDirectionZInput->text().toDouble());
    if (mWindow->isRunning()) mWindow->stopCalc();
    int particleIndex = mParticleIndexInput->text().toInt();
    mWindow->SetEnergyDefinitionAxis(particleIndex, direction, mStart, mEnd);
    QRect A = Bild->contentsRect();
    int w = A.width() - ScaleYWidth;
    double minR = xStart->text().toDouble(), maxR = xStop->text().toDouble();
    const Vector step = 0.01 * (mEnd - mStart), start = mStart + minR * step, end = mStart + maxR * step;
    mData.setParticleIndex(particleIndex);
    mData.Rescale(start, end, w + 1);
    if (mWindow->isRunning()) mWindow->stopCalc();
    mWindow->GetAxisEnergies(mData);
    Paint();
}

void PotentialDefiner::PSpektrum(QPainter &P, const QRect & A, bool /*PrintFN*/ )
{
    int w = A.width() - ScaleYWidth, h = A.height() - ScaleXHeight, di, n, N = mData.getNumnPoints();
    if (w <= 0 || h <= 0) return;
    double minR = xStart->text().toDouble(), maxR = xStop->text().toDouble();
    Vector start = mData.getEnd1(), end = mData.getEnd2();
    double Data[N], dataStart = (start - mStart).length(), length = (end - start).length(), step = length / (N-1), R, lastR;
    if (maxR < dataStart || minR > dataStart + length) return;
    int nStart = (minR > dataStart ? (minR - dataStart) / step : 0);
    step *= 100.0 / length;
    for (n=0; n<N; ++n) Data[n] = 0.0;
    for (n=0; n<7; n++)
    {
        if (n<6) P.setPen(CopyColor[n+1]);
        else P.setPen(CopyColor[0]);
        startLine(lastR = R = dataStart + nStart * step, Data[nStart]);
        for (di = nStart + 1; di < N && lastR < maxR; di++)
        {
            lastR = R;
            switch(n)
            {
            case 0:
                Data[di] += mData.GetFirstBound(di);
                break;
            case 1:
                Data[di] += mData.GetSecondBound(di);
                break;
            case 2:
                Data[di] += mData.GetThirdBound(di);
                break;
            case 3:
                Data[di] += mData.GetFourthBound(di);
                break;
            case 4:
                Data[di] += mData.GetFifthBound(di);
                break;
            case 5:
                Data[di] += mData.GetSixthBound(di);
                break;
            case 6:
                Data[di] += mData.GetUnbound(di);
                break;
            }
            continueLine(P, R += step, Data[di]);
        }
    }
    //DrawPoints(P, A);
}
