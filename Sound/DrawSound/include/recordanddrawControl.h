//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once

#include <QWidget>
#include <QAudioFormat>
#include <QAudioDecoder>


class DiagWindow;
class SoundMainWindow;
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
    enum Message{ChannelCountLargerOne, NotAllDataCouldBeWritten, DecodeError};

    struct AssignmentElement
    {
        double* data = nullptr;
        QString phoneme;
    };

    SoundRecordAndDrawControl(SoundMainWindow* MW);
    ~SoundRecordAndDrawControl();

    void Draw(const int sampleSize, const int sampleRate, const QAudioFormat::SampleType sampleType, const char* const inputData, const int nBytes);
    void Save(const char* const inputData, const int nBytes);
    void InitializeAssignmentData(const int numElements, const int numOscillators);
    AssignmentElement* GetAssignmentData(int& numElements) const;

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

    inline SoundMainWindow* GetMW() const
    {
        return mMW;
    }

private slots:
    void StartRecording();
    void Stop();
    void Draw();
    void Decode();
    void showInputFileDialog();
    void showOutputFileDialog();
    void SplitFileIntoPackets();
    void ReadyToDraw();
    void Error(QAudioDecoder::Error error);
    void BufferReady();
    void ShowMessage(Message message);

signals:
    void showMessage(Message);

private:
    enum DecodingFor{DF_Draw, DF_File, DF_Nothing};

    void VerifyFileExists(QString deviceName, QFile*& file, QLineEdit* edit);
    bool DetermineSampleTypeAndSize();
    void draw(const char* const inputData, const int nBytes);
    void writeRST();
    void createDecoder();
    void clearAssignmentData();

    DecodingFor mDecodingFor = DF_Nothing;
    QComboBox *mInputSelectorBox;
    QPushButton *mStartButton, *mStopButton, *mDrawButton, *mDecodeButton, *mInputFileDialogButton, *mOutputFileDialogButton, *mSplitFileButton;
    QLabel *mSizeDisplay, *mLengthDisplay;
    QLineEdit *mInputFileNameEdit, *mOutputFileNameEdit, *mPacketSizeEdit;
    QAudioInput* mInput;
    QAudioDecoder* mDecoder = nullptr;
    QByteArray mDecodeBuffer;
    const QString RST = "RST";
    QFile *mInputFile, *mOutputFile;
    QAudioFormat::SampleType mSampleType;
    int mSampleSize, mSampleRate, mNumChannels, mNumAssignmentElements;
    qint64 mProcessedUSec;
    std::vector<DiagWindow*> mFrequencyWindows;
    SoundMainWindow* mMW;
    AssignmentElement* mAssignmentElements = nullptr;
};
