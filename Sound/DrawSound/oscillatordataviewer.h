#pragma once

#include "DiagWindow.h"
#include "oscillatorarray.h"
#include "sounddrawwindow.h"


class SoundMainWindow;
class SoundWindow;
class QLineEdit;
class QComboBox;
class QKeyEvent;
class QPushButton;
class QElapsedTimer;


class OscillatorDataViewer : public DiagWindow
{
    Q_OBJECT
public:
    OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename, SoundWindow* soundWindow);
    virtual ~OscillatorDataViewer();

    void setLabels(const std::vector<SoundDrawWindow::Label>& labels);

private slots:
    void IncreaseTime();
    void DecreaseTime();

    void TimeChanged(const QString& value);
    void LabelChanged(int index);
    void KeyPressed(QKeyEvent *K);
    void ModeChanged(bool checked);

private:
    const OscillatorArray::Results& mData;
    std::vector<SoundDrawWindow::Label> mLabels;
    SoundWindow* mSoundWindow;
    QLineEdit* mTimeEdit, *mStepSizeEdit;
    QComboBox* mLabelBox;
    QPushButton* mModeButton;
    QElapsedTimer* mKeyTime;
    QString mKeyText;
    int mTimeIndex;
    double mHalfDeltaT;
};
