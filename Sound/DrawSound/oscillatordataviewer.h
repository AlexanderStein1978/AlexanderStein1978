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
    OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename);
    virtual ~OscillatorDataViewer();

    void setLabels(const std::vector<SoundDrawWindow::Label>& labels);

private slots:
    void IncreaseTime();
    void DecreaseTime();

    void TimeChanged(const QString& value);
    void LabelChanged(int index);

private:
    const OscillatorArray::Results& mData;
    std::vector<SoundDrawWindow::Label> mLabels;
    QLineEdit* mTimeEdit, *mStepSizeEdit;
    QComboBox* mLabelBox;
    int mTimeIndex;
    double mHalfDeltaT;
};
