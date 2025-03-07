#include "oscillatordiagram.h"
#include "soundmainwindow.h"
#include "oscillatorarray.h"
#include "utils.h"

#include <QLineEdit>


OscillatorDiagram::OscillatorDiagram(SoundMainWindow* MW) : DiagWindow(SimpleDiagWindow, MW), mTimePixelAssignments(nullptr), mFrequencyPixelAssignments(nullptr)
{
    setUnits("time [s]", "frequency [Hz]");
}

OscillatorDiagram::~OscillatorDiagram() noexcept
{
    if (nullptr != mTimePixelAssignments) delete[] mTimePixelAssignments;
    if (nullptr != mFrequencyPixelAssignments) delete[] mFrequencyPixelAssignments;
    if (nullptr != mData.data)
    {
        Destroy(mData.data, mData.numTimeSteps);
        delete[] mData.frequency;
        delete[] mData.time;
    }
}

void OscillatorDiagram::setData(OscillatorArray::Results& data)
{
    mData = data;
    mTimePixelAssignments = new int[mData.numTimeSteps];
    mFrequencyPixelAssignments = new int[OscillatorArray::NumOscillators];
	XMin = mData.time[0];
	XMax = mData.time[mData.numTimeSteps - 1];
	YMin = mData.frequency[0];
	YMax = mData.frequency[OscillatorArray::NumOscillators - 1];
	xStart->setText(QString::number(XMin, 'g', 11));
	xStop->setText(QString::number(XMax, 'g', 11));
	yStart->setText(QString::number(YMin));
	yStop->setText(QString::number(YMax));
    Paint();
}

void OscillatorDiagram::PSpektrum(QPainter& P, const QRect& A, bool PrintFN)
{
    if (XMax == XMin || YMax == YMin) return;
	QRgb *PP;
    double xmin = xStart->text().toDouble(), xmax = xStop->text().toDouble();
	double ymin = yStart->text().toDouble(), ymax = yStop->text().toDouble(), Tsc, Fsc, T;
    int BHeight = A.height() - ScaleTopOffset - ScaleXHeight;
	int BWidth = A.width() - ScaleYWidth, TIMin, TIMax, FIMin, FIMax, n, m, i, j;
    if (BHeight <= 0 || BWidth <= 0) return;
    for (TIMin = 0; TIMin < mData.numTimeSteps && mData.time[TIMin] < xmin; ++TIMin) mTimePixelAssignments[TIMin] = -1;
    if (TIMin == mData.numTimeSteps) return;
    for (TIMax = TIMin; TIMax < mData.numTimeSteps; ++TIMax) if (mData.time[TIMax] > xmax) break;
    if (TIMax == 0) return;
    int PictWidth = TIMax - TIMin;
    if (PictWidth > BWidth)
    {
        Tsc = static_cast<double>(PictWidth) / BWidth;
        for (n = TIMin, m=0, T = n; n < TIMax && m < BWidth; ++m, T += Tsc) while (n < mData.numTimeSteps && n <= T) mTimePixelAssignments[n++] = m;
        for (--m; n < TIMax; ++n) mTimePixelAssignments[n] = m;
        PictWidth = BWidth;
    }
    else for (n = TIMin; n < TIMax; ++n) mTimePixelAssignments[n] = n - TIMin;
    for (n = TIMax; n < mData.numTimeSteps; ++n) mTimePixelAssignments[n] = -1;
    for (FIMin = 0; FIMin < OscillatorArray::NumOscillators && mData.frequency[FIMin] < ymin; ++FIMin) mFrequencyPixelAssignments[FIMin] = -1;
    if (FIMin == OscillatorArray::NumOscillators) return;
    for (FIMax = FIMin; FIMax < OscillatorArray::NumOscillators; ++FIMax) if (mData.frequency[FIMax] > ymax) break;
    if (FIMax == 0) return;
    int PictHeight = FIMax - FIMin;
    if (PictHeight > BHeight)
    {
        Fsc = static_cast<double>(PictHeight) / BHeight;
        for (n = FIMin, m=0, T = n; n < FIMax && m < BHeight; ++m, T += Fsc) while (n < OscillatorArray::NumOscillators && n <= T)
            mFrequencyPixelAssignments[n++] = m;
        for (--m; n < FIMax; ++n) mFrequencyPixelAssignments[n] = m;
        PictHeight = BHeight;
    }
    else for (n = FIMin; n < FIMax; ++n) mFrequencyPixelAssignments[n] = n - FIMin;
    for (n = FIMax; n < OscillatorArray::NumOscillators; ++n) mFrequencyPixelAssignments[n] = -1;
    double **Pixel = Create(PictWidth, PictHeight);
	QImage *Pict = new QImage(PictWidth, PictHeight, QImage::Format_RGB32);
    double* MaxAmp = new double[PictHeight];
    for (n=0; n < PictHeight; ++n)
    {
        MaxAmp[n] = 0.0;
        for (m=0; m < PictWidth; ++m) Pixel[m][n] = 0.0;
    }
    for (n = TIMin; n < TIMax; ++n) for (m = FIMin; m < FIMax; ++m) if (mData.data[n][m] > Pixel[mTimePixelAssignments[n]][mFrequencyPixelAssignments[m]])
        Pixel[mTimePixelAssignments[n]][mFrequencyPixelAssignments[m]] = mData.data[n][m];
    for (n=0; n < PictWidth; ++n) for (m=0; m < PictHeight; ++m) if (Pixel[n][m] > MaxAmp[m]) MaxAmp[m] = Pixel[n][m];
    QRgb* PL;
    for (n=0, i = PictHeight - 1; n < PictHeight; ++n, --i)
    {
        double Isc = 1024.0 / MaxAmp[i];
        for (m=0, PL = reinterpret_cast<QRgb*>(Pict->scanLine(n)); m < PictWidth; ++m)
        {
            int value = static_cast<int>(Isc * Pixel[m][i]);
            if (value < 256) PL[m] = qRgb(0, 0, value);
            else if (value < 512)
            {
                value -= 256;
                PL[m] = qRgb(0, value, 255 - value);
            }
            else if (value < 768)
            {
                value -= 512;
                PL[m] = qRgb(value, 255 - value, 0);
            }
            else if (value < 1024)
            {
                value -= 768;
                PL[m] = qRgb(255, value, value);
            }
            else PL[m] = qRgb(255, 255, 255);
        }
    }
    Destroy(Pixel, PictWidth);
    delete[] MaxAmp;
	if (Image != nullptr) delete Image;
	Image = Pict;
	DiagWindow::PSpektrum(P, A, PrintFN);
}
