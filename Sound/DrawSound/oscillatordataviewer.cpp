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
    : DiagWindow(SimpleDiagWindow, MW, "Data files (*.dat)", ".dat", 1), mData(data), mTimeEdit(new QLineEdit("0.0", this))
    , mStepSizeEdit(new QLineEdit(QString::number(mData.time[1] - mData.time[0], 'g', 5))), mTimeIndex(0), mHalfDeltaT(0.5 * (data.time[1] - data.time[0]))
{
    setWindowTitle("Oscillator Data Viewer " + filename);
    setUnits("Frequency [Hz]", "Energy [arbitrary units]");
    QPushButton* increaseButton = new QPushButton("+", this), *decreaseButton = new QPushButton("-", this);
    SpektrumLayout->addWidget(new QLabel("Time [s]:", this), 0, 0);
    SpektrumLayout->addWidget(decreaseButton, 0, 1);
    SpektrumLayout->addWidget(mTimeEdit, 0, 2);
    SpektrumLayout->addWidget(increaseButton, 0, 3);
    SpektrumLayout->addWidget(new QLabel("step size:", this), 0, 4);
    SpektrumLayout->addWidget(mStepSizeEdit, 0, 5);
    mTimeEdit->setValidator(new QDoubleValidator(data.time[0], data.time[data.numTimeSteps - 1], 10, mTimeEdit));
    connect(increaseButton, SIGNAL(clicked()), this, SLOT(IncreaseTime()));
    connect(decreaseButton, SIGNAL(clicked()), this, SLOT(DecreaseTime()));
    connect(mTimeEdit, SIGNAL(textChanged(const QString&)), this, SLOT(TimeChanged(const QString&)));
}

void OscillatorDataViewer::DecreaseTime()
{
    mTimeEdit->setText(QString::number(mTimeEdit->text().toDouble() - mStepSizeEdit->text().toDouble(), 'f', 10));
}

void OscillatorDataViewer::IncreaseTime()
{
    mTimeEdit->setText(QString::number(mTimeEdit->text().toDouble() + mStepSizeEdit->text().toDouble(), 'f', 10));
}

void OscillatorDataViewer::TimeChanged(const QString& value)
{
    double newTime = value.toDouble();
    while (mTimeIndex < mData.numTimeSteps - 1 && mData.time[mTimeIndex] + mHalfDeltaT < newTime) ++mTimeIndex;
    while (mTimeIndex > 0 && mData.time[mTimeIndex] - mHalfDeltaT > newTime) --mTimeIndex;
    double** currentData = Create(OscillatorArray::NumOscillators, 2);
    for (int n=0; n < OscillatorArray::NumOscillators; ++n)
    {
        currentData[n][0] = mData.frequency[n];
        currentData[n][1] = mData.data[mTimeIndex][n];
    }
    setData(currentData, OscillatorArray::NumOscillators);
}
