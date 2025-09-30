//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "intensityhistogram.h"
#include "Spektrum.h"

#include <stdio.h>

#include <QPainter>
#include <QRect>

IntensityHistogram::IntensityHistogram(Spektrum* S, MainWindow* MW) : DiagWindow(IntensityHistogramPlot, MW)
{
    Spectrum = S;
    int n, N = 9 * width() / 10;
    S->getMinMaxIntensities(XMin, XMax);
    int *I = S->getIntensDist(N, XMin, XMax);
    for (YMin = YMax = 0.0, n=0; n<N; n++) if (I[n] > YMax) YMax = I[n];
    delete[] I;
    YMax *= 1.2;
    XMin = 0.0;
    XUnit = "Intensity";
    YUnit = "Number of peaks";
    xStart->setText("0.0");
    xStop->setText(QString::number(XMax, 'g', 6));
    yStart->setText("0");
    yStop->setText(QString::number(int(YMax)));
    xStartLabel->setText("Min intensity:");
    xStopLabel->setText("Max intensity:");
    yStartLabel->setText("Min count:");
    yStopLabel->setText("Max count:");
    setWindowTitle("Intensity histogram of spectrum " + S->getName());
}

void IntensityHistogram::PSpektrum(QPainter &P, const QRect &A, bool /*PrintFN*/)
{
    int left = ScaleYWidth, right = A.width(), x, n;
    int top = ScaleTopOffset, bottom = A.height() - ScaleXHeight, N = right - left + 1;
    double start = xStart->text().toDouble(), stop = xStop->text().toDouble(), R, M;
    int ystart = int(yStart->text().toDouble()), ystop = int(yStop->text().toDouble());
    int *I = Spectrum->getIntensDist(N, start, stop);
    for (x = left, n=0; n<N; n++, x++)
    {
        if (I[n] > ystart) P.drawLine(x, int(YO), x, (I[n] < ystop ? int(YO - YSF * double(I[n])) : top));
        //printf("I[%d]=%d, YO=%f, YSF=%f\n", n, I[n], YO, YSF);
    }
    Spectrum->getRMPH(R, M);
    if (R > start && R < stop)
    {
        P.setPen(QColor(255, 0, 0));
        x = int(XO + XSF * R);
        P.drawLine(x, top, x, bottom);
    }
    if (M > start && M < stop)
    {
        P.setPen(QColor(0, 0, 255));
        x = int(XO + XSF * M);
        P.drawLine(x, top, x, bottom);
    }
    delete[] I;
}
