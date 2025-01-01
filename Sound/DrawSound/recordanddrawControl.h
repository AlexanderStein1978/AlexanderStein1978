#pragma once

#include <QWidget>
#include <QAudioFormat>


class DiagWindow;
class QPushButton;
class QLabel;
class QComboBox;
class QAudioInput;
class QFile;
class QLineEdit;
class QAudioDecoder;


class SoundRecordAndDrawControl : public QWidget
{
    Q_OBJECT

public:
    SoundRecordAndDrawControl();
    ~SoundRecordAndDrawControl();

    inline void AddFrequencyWindow(DiagWindow* newWindow)
    {
        mFrequencyWindows.push_back(newWindow);
    }

    inline const std::vector<DiagWindow*>& GetTheFrequencyWindows() const
    {
        return mFrequencyWindows;
    }

    inline int GetNumberFrequencyWindows() const
    {
        return static_cast<int>(mFrequencyWindows.size());
    }

    inline DiagWindow* GetWindow(int index) const
    {
        return mFrequencyWindows[index];
    }

private slots:
    void StartRecording();
    void Stop();
    void Draw();
    void showFileDialog();
    void SplitFileIntoPackets();
    void ReadyToDraw();

private:
    void VerifyFileExists(QString deviceName);
    bool DetermineSampleTypeAndSize();
    void draw(const char* const inputData, const int nBytes);

    QComboBox *mInputSelectorBox;
    QPushButton *mStartButton, *mStopButton, *mDrawButton, *mFileDialogButton, *mSplitFileButton;
    QLabel* mSizeDisplay, *mLengthDisplay;
    QLineEdit *mFileNameEdit, *mPacketSizeEdit;
    QAudioInput* mInput;
    QAudioDecoder* mDecoder = nullptr;
    QFile* mFile;
    QAudioFormat::SampleType mSampleType;
    int mSampleSize, mSampleRate;
    qint64 mProcessedUSec;
    std::vector<DiagWindow*> mFrequencyWindows;
};
