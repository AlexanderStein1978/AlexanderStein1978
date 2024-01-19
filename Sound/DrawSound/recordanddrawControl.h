#pragma once

#include <QWidget>
#include <QAudioFormat>


class QPushButton;
class QLabel;
class QComboBox;
class QAudioInput;
class QFile;
class QLineEdit;


class SoundRecordAndDrawControl : public QWidget
{
    Q_OBJECT

public:
    SoundRecordAndDrawControl();
    ~SoundRecordAndDrawControl();

private slots:
    void StartRecording();
    void Stop();
    void Draw();
    void showFileDialog();

private:
    void VerifyFileExists(QString deviceName);
    bool DetermineSampleTypeAndSize();

    QComboBox *mInputSelectorBox;
    QPushButton *mStartButton, *mStopButton, *mDrawButton, *mFileDialogButton;
    QLabel* mSizeDisplay, *mLengthDisplay;
    QLineEdit* mFileNameEdit;
    QAudioInput* mInput;
    QFile* mFile;
    QAudioFormat::SampleType mSampleType;
    int mSampleSize, mSampleRate;
    qint64 mProcessedUSec;
};
