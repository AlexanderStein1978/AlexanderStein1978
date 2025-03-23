#pragma once

#include "DiagWindow.h"
#include "oscillatorarray.h"


class SoundMainWindow;
class QLineEdit;


class OscillatorDataViewer : public DiagWindow
{
    Q_OBJECT
public:
    OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename);

private slots:
    void IncreaseTime();
    void DecreaseTime();

    void TimeChanged(const QString& value);

private:
    const OscillatorArray::Results& mData;
    QLineEdit* mTimeEdit, *mStepSizeEdit;
    int mTimeIndex;
    double mHalfDeltaT;
};
