#include "recordanddrawControl.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QAudioDeviceInfo>
#include <QAudioInput>
#include <QBuffer>
#include <QMessageBox>

#include <algorithm>

#include "DiagWindow.h"
#include "utils.h"


namespace
{
    const QString SizeString("Size: "), LengthString("Length: ");
}


SoundRecordAndDrawControl::SoundRecordAndDrawControl() : mInputSelectorBox(new QComboBox(this)), mStartButton(new QPushButton("Start recording", this)),
    mStopAndDrawButton(new QPushButton("Stop recording and draw", this)), mSizeDisplay(new QLabel(SizeString, this)), mLengthDisplay(new QLabel(LengthString, this)), mInput(nullptr), mBuffer(nullptr)
{
    setWindowTitle("Sound Record and Draw Control");
    QGridLayout *L = new QGridLayout(this);
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    for (QAudioDeviceInfo info : deviceList) mInputSelectorBox->addItem(info.deviceName());
    L->setColumnStretch(0, 1);
    L->setColumnStretch(1, 1);
    L->addWidget(new QLabel("Input device:", this), 0, 0);
    L->addWidget(mInputSelectorBox, 0, 1);
    L->addWidget(mStartButton, 1, 0);
    L->addWidget(mStopAndDrawButton, 1, 1);
    mStopAndDrawButton->setEnabled(false);
    L->addWidget(mSizeDisplay, 2, 0);
    L->addWidget(mLengthDisplay, 2, 1);
    connect(mStartButton, SIGNAL(clicked()), this, SLOT(StartRecording()));
    connect(mStopAndDrawButton, SIGNAL(clicked()), this, SLOT(StopAndDraw()));
}

SoundRecordAndDrawControl::~SoundRecordAndDrawControl()
{
    if (nullptr != mInput) delete mInput;
    if (nullptr != mBuffer) delete mBuffer;
}

void SoundRecordAndDrawControl::StartRecording()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    int deviceIndex = mInputSelectorBox->currentIndex();
    if (nullptr == mBuffer)
    {
        mBuffer = new QBuffer;
    }
    else
    {
        mBuffer->close();
        mBuffer->buffer().clear();
    }
    mBuffer->buffer().reserve(100000000);
    mBuffer->open(QIODevice::WriteOnly);
    QAudioFormat format;
    int sampleRate = 0;
    QList<int> ssrs = deviceList[deviceIndex].supportedSampleRates();
    for (int s : ssrs) if (s > sampleRate) sampleRate = s;
    format.setSampleRate(sampleRate);
    format.setChannelCount(1);
    mSampleSize = 0;
    QList<int> ssss = deviceList[deviceIndex].supportedSampleSizes();
    for (int s : ssss) if (s > mSampleSize) mSampleSize = s;
    if (0 == mSampleSize)
    {
        QMessageBox::warning(this, "DrawSound", "For the selected input device the needed information could not be found!");
        return;
    }
    format.setSampleSize(mSampleSize);
    format.setCodec("audio/pcm");
    format.setByteOrder(QAudioFormat::LittleEndian);
    mSampleType = QAudioFormat::Unknown;
    QList<QAudioFormat::SampleType> ssts = deviceList[deviceIndex].supportedSampleTypes();
    for (QAudioFormat::SampleType t : ssts)
    {
        switch (mSampleType)
        {
            case QAudioFormat::Unknown:
               mSampleType = t;
               continue;
            case QAudioFormat::SignedInt:
                if (t == QAudioFormat::Float) mSampleType = t;
                continue;
            case QAudioFormat::UnSignedInt:
                if (t == QAudioFormat::SignedInt || t == QAudioFormat::Float) mSampleType = t;
                continue;
            case QAudioFormat::Float:
                // do nothing
                break;
            default:
                mSampleType = t;
                continue;
                break;
        }
        break;
    }
    format.setSampleType(mSampleType);
    if (nullptr != mInput) delete mInput;
    mInput = new QAudioInput(deviceList[deviceIndex], format, this);
    mInput->start(mBuffer);
    mInputSelectorBox->setEnabled(false);
    mStartButton->setEnabled(false);
    mStopAndDrawButton->setEnabled(true);
}

void SoundRecordAndDrawControl::StopAndDraw()
{
    mInput->suspend();
    int nBytes = mInput->bytesReady();
    int nSamples = nBytes / mSampleSize * 8, i, pos;
    if (0 < nSamples)
    {
        double passedTime = 0.000001 * mInput->processedUSecs(), currentTime;
        double timeStep = passedTime / nSamples;
        QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
        int deviceIndex = mInputSelectorBox->currentIndex();
        double **data = Create(nSamples, 2);
        QByteArray inputArray = mBuffer->data();
        char* inputData = inputArray.data();
        for (i=pos=0, currentTime = 0.0; i < nSamples; ++i, currentTime += timeStep, pos += (mSampleSize / 8))
        {
            data[i][0] = currentTime;
            switch (mSampleType)
            {
                case QAudioFormat::UnSignedInt:
                    switch (mSampleSize)
                    {
                        case 8:
                            data[i][1] = static_cast<u_int8_t>(inputData[i]);
                            break;
                        case 16:
                            data[i][1] = *reinterpret_cast<u_int16_t*>(inputData + pos);
                            break;
                        default:
                            u_int32_t sample = 0;
                            memcpy(&sample, inputData + pos, mSampleSize / 8);
                            data[i][1] = sample;
                            break;
                    }
                case QAudioFormat::Float:
                    if (mSampleSize < 32)
                    {
                        _Float32 sample = 0.0f;
                        memcpy(&sample, inputData + pos, mSampleSize / 8);
                        data[i][1] = sample;
                    }
                    else if (mSampleSize == 32) data[i][1] = *reinterpret_cast<_Float32*>(inputData + pos);
                    else if (mSampleSize < 64)
                    {
                        double sample = 0.0;
                        memcpy(&sample, inputData + pos, mSampleSize / 8);
                        data[i][1] = sample;
                    }
                    else if (mSampleSize == 64) data[i][1] = *reinterpret_cast<double*>(inputData + pos);
                    break;
                case QAudioFormat::Unknown:
                case QAudioFormat::SignedInt:
                default:
                    switch (mSampleSize)
                    {
                        case 8:
                            data[i][1] = inputData[i];
                            break;
                        case 16:
                            data[i][1] = *reinterpret_cast<int16_t*>(inputData + i);
                            break;
                        default:
                            int32_t sample = 0;
                            memcpy(&sample, inputData + i, mSampleSize);
                            data[i][1] = sample;
                            break;
                    }
            }
        }
        DiagWindow* window = new DiagWindow;
        window->setData(data, nSamples);
        window->setUnits("time [s]", "intensity");
        window->show();
    }
    else QMessageBox::information(this, "DrawSound", "With the selected input device no sound could be recorded!");
    mInput->stop();
    mInputSelectorBox->setEnabled(true);
    mStartButton->setEnabled(true);
    mStopAndDrawButton->setEnabled(false);
}
