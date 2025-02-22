#include "oscillatordiagram.h"
#include "soundmainwindow.h"
#include "oscillatorarray.h"


OscillatorDiagram::OscillatorDiagram(SoundMainWindow* MW) : DiagWindow(MW), mTimePixelAssignments(nullptr), mFrequencyPixelAssignments(nullptr)
{
    setUnits("time [s]", "frequency [Hz]");
}

OscillatorDiagram::~OscillatorDiagram() noexcept
{
    if (nullptr != mTimePixelAssignments) delete[] mTimePixelAssignments;
    if (nullptr != mFrequencyPixelAssignments) delete[] mFrequencyPixelAssignments;
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
	double ymin = yStart->text().toDouble(), ymax = yStop->text().toDouble(), Tsc, Fsc, MaxAmp = 0.0, T;
    int BHeight = A.height() - ScaleTopOffset - ScaleXHeight;
	int BWidth = A.width() - ScaleYWidth, TIMin, TIMax, FIMin, FIMax, n, m, i, j;
    for (TIMin = 0; TIMin < mData.numTimeSteps && mData.time[TIMin] < xmin; ++TIMin) mTimePixelAssignments[TIMin] = -1;
    if (TIMin == mData.numTimeSteps) return;
    for (TIMax = TIMin; TIMax < mData.numTimeSteps; ++TIMax) if (mData.time[TIMax] > xmax) break;
    if (TIMax == 0) return;
    int PictWidth = TIMax - TIMin + 1;
    if (PictWidth > BWidth)
    {
        Tsc = PictWidth / BWidth;
        for (n = TIMin, m=0, T = mData.time[n]; n <= TIMax; ++m, T += Tsc) while (n < mData.numTimeSteps && mData.time[n] <= T) mTimePixelAssignments[n++] = m;
        PictWidth = BWidth;
    }
    else
    {
        for (n = TIMin; n <= TIMax; ++n) mTimePixelAssignments[n] = n - TIMin;
        Tsc = 1.0;
    }
    for (n = TIMax + 1; n < mData.numTimeSteps; ++n) mTimePixelAssignments[n] = -1;
    for (FIMin = 0; FIMin < OscillatorArray::NumOscillators && mData.frequency[FIMin] < ymin; ++FIMin) mFrequencyPixelAssignments[FIMin] = -1;
    if (FIMin == OscillatorArray::NumOscillators) return;
    for (FIMax = FIMin; FIMax < OscillatorArray::NumOscillators; ++FIMax) if (mData.frequency[FIMax] > ymax) break;
    if (FIMax == 0) return;
    int PictHeight = FIMax - FIMin + 1;
    if (PictHeight > BHeight)
    {
        Fsc = PictHeight / BHeight;
        for (n = FIMin, m=0, T = mData.frequency[n]; n <= FIMax; ++m, T += Tsc) while (n < OscillatorArray::NumOscillators && mData.frequency[n] <= T)
            mFrequencyPixelAssignments[n++] = m;
        PictHeight = BHeight;
    }
    else
    {
        for (n = FIMin; n <= FIMax; ++n) mFrequencyPixelAssignments[n] = n - FIMin;
        Fsc = 1.0;
    }
    for (n = FIMax + 1; n < OscillatorArray::NumOscillators; ++n) mFrequencyPixelAssignments[n] = -1;
    double **Pixel = Create(PictWidth, PictHeight);
	QImage *Pict = new QImage(PictWidth, PictHeight, QImage::Format_RGB32);
    for (n=0, i = TIMin; n < PictWidth; ++n) for (m=0, j = FIMin; m < PictHeight; ++m)
    {
        for (Pixel[n][m] = 0.0; i < mData.numTimeSteps && mTimePixelAssignments[i] == n; ++i)
            for ( ; j < OscillatorArray::NumOscillators && mFrequencyPixelAssignments[j] == m) if (mData.data[i][j] > Pixel[n][m])
                Pixel[n][m] = mData.data[i][j];
        if (Pixel[n][m] > MaxAmp) MaxAmp = Pixel[n][m];
    }

}
