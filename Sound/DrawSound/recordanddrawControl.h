#pragma once

#include <QWidget>
#include <QAudioFormat>


class QPushButton;
class QLabel;
class QComboBox;
class QAudioInput;
class QBuffer;


class SoundRecordAndDrawControl : public QWidget
{
    Q_OBJECT

public:
    SoundRecordAndDrawControl();
    ~SoundRecordAndDrawControl();

private slots:
    void StartRecording();
    void StopAndDraw();

private:
    QComboBox *mInputSelectorBox;
    QPushButton *mStartButton, *mStopAndDrawButton;
    QLabel* mSizeDisplay, *mLengthDisplay;
    QAudioInput* mInput;
    QBuffer* mBuffer;
    QAudioFormat::SampleType mSampleType;
    int mSampleSize;
};
