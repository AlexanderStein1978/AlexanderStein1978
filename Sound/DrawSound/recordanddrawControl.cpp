#include "recordanddrawControl.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QAudioDeviceInfo>
#include <QAudioInput>
#include <QFile>
#include <QFileDialog>
#include <QDataStream>
#include <QMessageBox>

#include <algorithm>

#include "soundwindow.h"
#include "utils.h"


namespace
{
    const QString SizeString("Size: "), LengthString("Length: "), DefaultDirectory(DATA_DIRECTORY "/Recordings/");
}


SoundRecordAndDrawControl::SoundRecordAndDrawControl() : mInputSelectorBox(new QComboBox(this)), mStartButton(new QPushButton("Start recording", this)),
    mStopButton(new QPushButton("Stop recording", this)), mDrawButton(new QPushButton("Draw", this)), mFileDialogButton(new QPushButton("...", this)), mSizeDisplay(new QLabel(SizeString, this)),
    mLengthDisplay(new QLabel(LengthString, this)), mFileNameEdit(new QLineEdit(this)), mInput(nullptr), mFile(nullptr), mSampleType(QAudioFormat::Unknown), mSampleSize(0), mSampleRate(0),
    mProcessedUSec(0u)
{
    setWindowTitle("Sound Record and Draw Control");
    QGridLayout *L = new QGridLayout(this);
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    for (QAudioDeviceInfo info : deviceList) mInputSelectorBox->addItem(info.deviceName());
    L->addWidget(new QLabel("Input device:", this), 0, 0);
    L->addWidget(mInputSelectorBox, 0, 1, 1, 2);
    L->addWidget(new QLabel("Filename:", this), 1, 0);
    QGridLayout *FL  = new QGridLayout;
    L->addLayout(FL, 1, 1, 1, 2);
    FL->addWidget(mFileNameEdit, 0, 0);
    FL->addWidget(mFileDialogButton, 0, 1);
    L->addWidget(mStartButton, 2, 0);
    L->addWidget(mStopButton, 2, 1);
    L->addWidget(mDrawButton, 2, 2);
    mStopButton->setEnabled(false);
    L->addWidget(mSizeDisplay, 3, 0);
    L->addWidget(mLengthDisplay, 3, 1);
    connect(mStartButton, SIGNAL(clicked()), this, SLOT(StartRecording()));
    connect(mStopButton, SIGNAL(clicked()), this, SLOT(Stop()));
    connect(mDrawButton, SIGNAL(clicked()), this, SLOT(Draw()));
    connect(mFileDialogButton, SIGNAL(clicked()), this, SLOT(showFileDialog()));
}

SoundRecordAndDrawControl::~SoundRecordAndDrawControl()
{
    if (nullptr != mInput) delete mInput;
    if (nullptr != mFile) delete mFile;
}

void SoundRecordAndDrawControl::showFileDialog()
{
    QString fileName = mFileNameEdit->text();
    if (!fileName.contains(QRegExp("[/\\\\]"))) fileName = DefaultDirectory + fileName;
    fileName = QFileDialog::getSaveFileName(this, "Select filename", fileName);
    if (!fileName.isEmpty()) mFileNameEdit->setText(fileName);
}

void SoundRecordAndDrawControl::VerifyFileExists(QString deviceName)
{
    QString filename = mFileNameEdit->text();
    if (filename.isEmpty()) filename = deviceName.replace(QRegularExpression("[.,:=]"), "_");
    if (!filename.contains('.')) filename += ".dat";
    if (!filename.contains(QRegExp("[/\\\\]"))) filename = DefaultDirectory + filename;
    if (nullptr != mFile) delete mFile;
    mFile = new QFile(filename, this);
}

bool SoundRecordAndDrawControl::DetermineSampleTypeAndSize()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    int deviceIndex = mInputSelectorBox->currentIndex();
    mSampleSize = 0;
    QList<int> ssss = deviceList[deviceIndex].supportedSampleSizes();
    for (int s : ssss) if (s > mSampleSize) mSampleSize = s;
    if (0 == mSampleSize)
    {
        QMessageBox::warning(this, "DrawSound", "For the selected input device the needed information could not be found!");
        return false;
    }
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
    mSampleRate = 0;
    QList<int> ssrs = deviceList[deviceIndex].supportedSampleRates();
    for (int s : ssrs) if (s > mSampleRate) mSampleRate = s;
    mProcessedUSec = 0u;
    return true;
}

void SoundRecordAndDrawControl::StartRecording()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    int deviceIndex = mInputSelectorBox->currentIndex();
    VerifyFileExists(mInputSelectorBox->currentText());
    mFile->open(QIODevice::WriteOnly);
    QAudioFormat format;
    if (!DetermineSampleTypeAndSize()) return;
    format.setSampleRate(mSampleRate);
    format.setChannelCount(1);
    format.setSampleSize(mSampleSize);
    format.setCodec("audio/pcm");
    format.setByteOrder(QAudioFormat::LittleEndian);
    format.setSampleType(mSampleType);
    if (nullptr != mInput) delete mInput;
    mInput = new QAudioInput(deviceList[deviceIndex], format, this);
    mInput->start(mFile);
    mInputSelectorBox->setEnabled(false);
    mStartButton->setEnabled(false);
    mStopButton->setEnabled(true);
}

void SoundRecordAndDrawControl::Stop()
{
    mProcessedUSec = mInput->processedUSecs();
    mInput->stop();
    mFile->close();
    mInputSelectorBox->setEnabled(true);
    mStartButton->setEnabled(true);
    mStopButton->setEnabled(false);
}

void SoundRecordAndDrawControl::Draw()
{
    VerifyFileExists(mInputSelectorBox->currentText());
    int nBytes = mFile->size();
    mFile->open(QIODevice::ReadOnly);
    QDataStream stream(mFile);
    char* inputData = new char[nBytes];
    nBytes = stream.readRawData(inputData, nBytes);
    if (0 == mSampleSize && !DetermineSampleTypeAndSize()) return;
    int nSamples = nBytes / mSampleSize * 8, i, pos;
    if (nullptr != inputData && 0 < nBytes)
    {
        if (0u == mProcessedUSec) mProcessedUSec = static_cast<double>(nBytes) * 8000000 / (mSampleSize * mSampleRate);
        double passedTime = 0.000001 * mProcessedUSec, currentTime;
        double timeStep = passedTime / nSamples;
        double **data = Create(nSamples, 2);
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
                    // printf("(%g|%g), ", data[i][0], data[i][1]);
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
                            data[i][1] = *reinterpret_cast<int16_t*>(inputData + pos);
                            break;
                        default:
                            int32_t sample = 0;
                            memcpy(&sample, inputData + pos, mSampleSize / 8);
                            data[i][1] = sample;
                            break;
                    }
            }
        }
        SoundDrawWindow* window = new SoundWindow(this, mFileNameEdit->text(), mSampleRate);
        window->setData(data, nSamples);
        window->show();
    }
    else QMessageBox::information(this, "DrawSound", "With the selected input device no sound could be recorded!");
    if (nullptr != inputData) delete[] inputData;
}
