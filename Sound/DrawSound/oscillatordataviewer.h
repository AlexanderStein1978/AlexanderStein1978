#pragma once

#include "DiagWindow.h"
#include "oscillatorarray.h"
#include "sounddrawwindow.h"


class SoundMainWindow;
class QLineEdit;
class QComboBox;


class OscillatorDataViewer : public DiagWindow
{
    Q_OBJECT
public:
    OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename, const std::vector<SoundDrawWindow::Label>& labels);
    ~OscillatorDataViewer();

private slots:
    void IncreaseTime();
    void DecreaseTime();

    void TimeChanged(const QString& value);
    void LabelChanged(int index);

private:
    const OscillatorArray::Results& mData;
    const std::vector<SoundDrawWindow::Label>& mLabels;
    int* mLabelIndices;
    QLineEdit* mTimeEdit, *mStepSizeEdit;
    QComboBox* mLabelBox;
    int mTimeIndex;
    double mHalfDeltaT;
};
