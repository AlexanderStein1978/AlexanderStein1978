//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once

#include <QWidget>
#include <QAudioFormat>
#include <QAudioDecoder>
<<<<<<< HEAD
#include <QMediaCaptureSession>
=======
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877


class DiagWindow;
class SoundMainWindow;
class QPushButton;
class QLabel;
class QComboBox;
class QAudioInput;
class QFile;
class QLineEdit;
<<<<<<< HEAD
class QMediaRecorder;
class QMediaDevices;
=======
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877


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

<<<<<<< HEAD
    void Draw(const int sampleSize, const int sampleRate, const QAudioFormat::SampleFormat sampleType, const char* const inputData, const int nBytes);
=======
    void Draw(const int sampleSize, const int sampleRate, const QAudioFormat::SampleType sampleType, const char* const inputData, const int nBytes);
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
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
<<<<<<< HEAD
	void updateFormats();
=======
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

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

<<<<<<< HEAD
	bool mUpdatingFormats = false;
    DecodingFor mDecodingFor = DF_Nothing;
    QComboBox *mInputSelectorBox, *mFileFormatBox, *mCodecBox;
=======
    DecodingFor mDecodingFor = DF_Nothing;
    QComboBox *mInputSelectorBox;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    QPushButton *mStartButton, *mStopButton, *mDrawButton, *mDecodeButton, *mInputFileDialogButton, *mOutputFileDialogButton, *mSplitFileButton;
    QLabel *mSizeDisplay, *mLengthDisplay;
    QLineEdit *mInputFileNameEdit, *mOutputFileNameEdit, *mPacketSizeEdit;
    QAudioInput* mInput;
    QAudioDecoder* mDecoder = nullptr;
    QByteArray mDecodeBuffer;
    const QString RST = "RST";
    QFile *mInputFile, *mOutputFile;
<<<<<<< HEAD
    QAudioFormat::SampleFormat mSampleType;
=======
    QAudioFormat::SampleType mSampleType;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    int mSampleSize, mSampleRate, mNumChannels, mNumAssignmentElements;
    qint64 mProcessedUSec;
    std::vector<DiagWindow*> mFrequencyWindows;
    SoundMainWindow* mMW;
    AssignmentElement* mAssignmentElements = nullptr;
<<<<<<< HEAD
	QMediaCaptureSession m_captureSession;
    QMediaRecorder *m_audioRecorder = nullptr;
    QMediaDevices *m_mediaDevices = nullptr;
=======
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
};
