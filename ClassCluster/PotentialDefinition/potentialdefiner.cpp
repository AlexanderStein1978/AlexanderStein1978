#include "potentialdefiner.h"
#include "window.h"

#include <QLineEdit>
#include <QPainter>
#include <QGridLayout>
#include <QLabel>


PotentialDefiner::PotentialDefiner(Window* window, MainWindow* MW) : DiagWindow(SimpleDiagWindow, MW, "", "*.*", 1), mWindow(window), mParticleIndexInput(new QLineEdit("0", this)),
    mDirectionXInput(new QLineEdit("0.0", this)), mDirectionYInput(new QLineEdit("0.0", this)), mDirectionZInput(new QLineEdit("0.0", this))
{
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
    SpektrumLayout->addLayout(layout, 0, 0, 1, 9);
    connect(mParticleIndexInput, SIGNAL(editingFinished()), this, SLOT(InputChanged()));
    connect(mDirectionXInput, SIGNAL(editingFinished()), this, SLOT(InputChanged()));
    connect(mDirectionYInput, SIGNAL(editingFinished()), this, SLOT(InputChanged()));
    connect(mDirectionZInput, SIGNAL(editingFinished()), this, SLOT(InputChanged()));
}

void PotentialDefiner::PSpektrum(QPainter &P, const QRect & A, bool /*PrintFN*/ )
{
    int i, w = A.width() - ScaleYWidth, h = A.height() - ScaleXHeight, di;
    if (w <= 0 || h <= 0) return;
    int l = A.left() + ScaleYWidth - 1;
    int /*t = A.top(), r = l + w,*/ lp, n;
    int p=0, b = A.bottom() - ScaleXHeight;
    double minR = xStart->text().toDouble(), maxR = xStop->text().toDouble();
    double minE = yStart->text().toDouble(), maxE = yStop->text().toDouble();
    double ESc = (double)h / (maxE - minE), Data[w];
    const Vector diff = mStart - mEnd, start = mStart + minR * diff, end = mStart + maxR * diff;
    mData.Rescale(start, end, w + 1);
    mWindow->GetAxisEnergies(mData);
    for (i=0; i<w; ++i) Data[i] = 0.0;
    for (n=0; n<6; n++)
    {
        P.setPen(CopyColor[n]);
        for (di=0, i=l; i <= l+w && di <= w; i++, di++)
        {
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
                Data[di] += mData.GetSecondOrderBound(di);
                break;
            case 5:
                Data[di] += mData.GetUnbound(di);
                break;
            }
            lp = p;
            p = b - int(ESc * (Data[di] - minE)) + 1;
            if (p > b) p = b + 1;
            if (i>l ? (p > lp ? p - lp : lp - p) > 1 : false) P.drawLine(i-1, lp, i, p);
            else P.drawPoint(i, p);
        }
    }
    //DrawPoints(P, A);
}
