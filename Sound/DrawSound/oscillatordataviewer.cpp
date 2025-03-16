#include "oscillatordataviewer.h"
#include "soundmainwindow.h"
#include "utils.h"

#include <QLineEdit>
#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QDoubleValidator>

#include <cmath>


OscillatorDataViewer::OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename)
    : DiagWindow(SimpleDiagWindow, MW, "Data files (*.dat)", ".dat", 1), mData(data), mTimeEdit(new QLineEdit("0.0", this)), mTimeIndex(0)
    , mHalfDeltaT(0.5 * (data.time[1] - data.time[0]))
{
    setWindowTitle("Oscillator Data Viewer " + filename);
    setUnits("Frequency [Hz]", "Energy [arbitrary units]");
    QPushButton* increaseButton = new QPushButton("+", this), *decreaseButton = new QPushButton("-", this);
    SpektrumLayout->addWidget(new QLabel("Time [s]:", this), 0, 0);
    SpektrumLayout->addWidget(decreaseButton, 0, 1);
    SpektrumLayout->addWidget(mTimeEdit, 0, 2);
    SpektrumLayout->addWidget(increaseButton, 0, 3);
    mTimeEdit->setValidator(new QDoubleValidator(data.time[0], data.time[data.numTimeSteps - 1], 10, mTimeEdit));
    connect(increaseButton, SIGNAL(clicked()), this, SLOT(IncreaseTime()));
    connect(decreaseButton, SIGNAL(clicked()), this, SLOT(DecreaseTime()));
    connect(mTimeEdit, SIGNAL(textChanged(const QString&)), this, SLOT(TimeChanged(const QString&)));
}

void OscillatorDataViewer::DecreaseTime()
{
    if (mTimeIndex > 0) mTimeEdit->setText(QString::number(mData.time[--mTimeIndex], 'f', 10));
}

void OscillatorDataViewer::IncreaseTime()
{
    if (mTimeIndex < mData.numTimeSteps - 1) mTimeEdit->setText(QString::number(mData.time[++mTimeIndex], 'f', 10));
}

void OscillatorDataViewer::TimeChanged(const QString& value)
{
    double newTime = value.toDouble();
    while (mData.time[mTimeIndex] + mHalfDeltaT < newTime) ++mTimeIndex;
    while (mData.time[mTimeIndex] - mHalfDeltaT > newTime) --mTimeIndex;
    double** currentData = Create(OscillatorArray::NumOscillators, 2);
    for (int n=0; n < OscillatorArray::NumOscillators; ++n)
    {
        currentData[n][0] = mData.frequency[n];
        currentData[n][1] = mData.data[mTimeIndex][n];
    }
    setData(currentData, OscillatorArray::NumOscillators);
}
