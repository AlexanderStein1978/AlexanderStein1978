//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "recordanddrawControl.h"

#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QAudioDevice>
#include <QMediaDevices>
#include <QAudioInput>
#include <QAudioDecoder>
#include <QFile>
#include <QFileDialog>
#include <QDataStream>
#include <QMessageBox>
#include <QtEndian>
#include <QMediaRecorder>
#include <QMediaFormat>

#include "soundwindow.h"
#include "soundmainwindow.h"
#include "utils.h"
#include "definesamplesizedialog.h"


namespace
{
    const QString SizeString("Size: "), LengthString("Length: "), DefaultDirectory(DATA_DIRECTORY "/Recordings/");
}


SoundRecordAndDrawControl::SoundRecordAndDrawControl(SoundMainWindow* MW) : mInputSelectorBox(new QComboBox(this)), mFileFormatBox(new QComboBox(this)), mCodecBox(new QComboBox(this)),
	mStartButton(new QPushButton("Start recording", this)), mStopButton(new QPushButton("Stop recording", this)), mDrawButton(new QPushButton("Draw", this)), mDecodeButton(new QPushButton("Decode", this)),
	mInputFileDialogButton(new QPushButton("...", this)), mOutputFileDialogButton(new QPushButton("...", this)), mSplitFileButton(new QPushButton("Split recording", this)), mSizeDisplay(new QLabel(SizeString, this)),
    mLengthDisplay(new QLabel(LengthString, this)), mInputFileNameEdit(new QLineEdit(this)), mOutputFileNameEdit(new QLineEdit(this)), mPacketSizeEdit(new QLineEdit("10", this)), mInput(nullptr), mInputFile(nullptr), mOutputFile(nullptr), mSampleType(QAudioFormat::Unknown), mSampleSize(0), mSampleRate(0), mProcessedUSec(0u), mMW(MW)
{
    setWindowTitle("Sound Record and Draw Control");
    QGridLayout *L = new QGridLayout(this);
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    for (QAudioDevice device : deviceList) mInputSelectorBox->addItem(device.description(), QVariant::fromValue(device));
    L->addWidget(new QLabel("Input device:", this), 0, 0);
    L->addWidget(mInputSelectorBox, 0, 1, 1, 3);
	mInputSelectorBox->setEditable(false);
	L->addWidget(new QLabel("File format:", this), 1, 0);
	L->addWidget(mFileFormatBox, 1, 1, 1, 3);
	mFileFormatBox->setEditable(false);
	L->addWidget(new QLabel("Codec:", this), 2, 0);
	L->addWidget(mCodecBox, 2, 1, 1, 3);
	mCodecBox->setEditable(false);
    L->addWidget(new QLabel("Input filename:", this), 3, 0);
    L->addWidget(new QLabel("Output filename:", this), 4, 0);
    QGridLayout *FL  = new QGridLayout;
    L->addLayout(FL, 3, 1, 2, 2);
    FL->addWidget(mInputFileNameEdit, 0, 0);
    FL->addWidget(mInputFileDialogButton, 0, 1);
    FL->addWidget(mOutputFileNameEdit, 1, 0);
    FL->addWidget(mOutputFileDialogButton, 1, 1);
    L->addWidget(mStartButton, 5, 0);
    L->addWidget(mStopButton, 5, 1);
    L->addWidget(mDrawButton, 5, 2);
    L->addWidget(mDecodeButton, 5, 3);
    mStopButton->setEnabled(false);
    L->addWidget(mSizeDisplay, 6, 0);
    L->addWidget(mLengthDisplay, 6, 1);
    L->addWidget(mSplitFileButton, 7, 0);
    L->addWidget(new QLabel("Size of split files", this),7, 1);
    L->addWidget(mPacketSizeEdit, 7, 2, 1, 2);
    connect(mStartButton, SIGNAL(clicked()), this, SLOT(StartRecording()));
    connect(mStopButton, SIGNAL(clicked()), this, SLOT(Stop()));
    connect(mDrawButton, SIGNAL(clicked()), this, SLOT(Draw()));
    connect(mDecodeButton, SIGNAL(clicked()), this, SLOT(Decode()));
    connect(mInputFileDialogButton, SIGNAL(clicked()), this, SLOT(showInputFileDialog()));
    connect(mOutputFileDialogButton, SIGNAL(clicked()), this, SLOT(showOutputFileDialog()));
    connect(mSplitFileButton, SIGNAL(clicked()), this, SLOT(SplitFileIntoPackets()));
    connect(this, SIGNAL(showMessage(Message)), this, SLOT(ShowMessage(Message)), Qt::QueuedConnection);
	connect(mInputSelectorBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFormats()));
	connect(mFileFormatBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFormats()));
    m_audioRecorder = new QMediaRecorder(this);
    m_captureSession.setRecorder(m_audioRecorder);
    m_captureSession.setAudioInput(new QAudioInput(this));
	updateFormats();
}

SoundRecordAndDrawControl::~SoundRecordAndDrawControl()
{
    if (nullptr != mInput) delete mInput;
    if (nullptr != mInputFile) delete mInputFile;
    if (nullptr != mOutputFile) delete mOutputFile;
    if (nullptr != mDecoder)
    {
        disconnect(mDecoder);
        delete mDecoder;
    }
    clearAssignmentData();
}

void SoundRecordAndDrawControl::updateFormats()
{
    if (mUpdatingFormats) return;
    mUpdatingFormats = true;

    QMediaFormat format;
    if (mFileFormatBox->count()) format.setFileFormat(mFileFormatBox->currentData().value<QMediaFormat::FileFormat>());

    int currentIndex = 0;
    mFileFormatBox->clear();

    const QList<QMediaFormat::FileFormat> supportedFormats = format.supportedFileFormats(QMediaFormat::Encode);
    for (auto container : supportedFormats)
	{
        if (container < QMediaFormat::Mpeg4Audio) continue;
        if (container == format.fileFormat()) currentIndex = mFileFormatBox->count();
        mFileFormatBox->addItem(QMediaFormat::fileFormatDescription(container), QVariant::fromValue(container));
    }
    mFileFormatBox->setCurrentIndex(currentIndex);

    QMediaFormat::AudioCodec currentCodec = QMediaFormat::AudioCodec::Unspecified;
    if (mCodecBox->count()) currentCodec = mCodecBox->currentData().value<QMediaFormat::AudioCodec>();

    currentIndex = 0;
    mCodecBox->clear();
 
	const QList<QMediaFormat::AudioCodec> supportedCodecs = format.supportedAudioCodecs(QMediaFormat::Encode);
    for (auto codec : supportedCodecs)
	{
        if (codec == currentCodec) currentIndex = mCodecBox->count();
        mCodecBox->addItem(QMediaFormat::audioCodecDescription(codec), QVariant::fromValue(codec));
    }
    mCodecBox->setCurrentIndex(currentIndex);

    mUpdatingFormats = false;
}

void SoundRecordAndDrawControl::showInputFileDialog()
{
    QString filename = mInputFileNameEdit->text();
    if (!filename.contains(QRegularExpression("[/\\\\]"))) filename = DefaultDirectory + filename;
    filename = QFileDialog::getOpenFileName(this, "Select filename to write to", filename);
    if (!filename.isEmpty()) mInputFileNameEdit->setText(filename);
}

void SoundRecordAndDrawControl::showOutputFileDialog()
{
    QString fileName = mOutputFileNameEdit->text();
    if (!fileName.contains(QRegularExpression("[/\\\\]"))) fileName = DefaultDirectory + fileName;
    fileName = QFileDialog::getSaveFileName(this, "Select filename to open", fileName);
    if (!fileName.isEmpty()) mOutputFileNameEdit->setText(fileName);
}

void SoundRecordAndDrawControl::VerifyFileExists(QString deviceName, QFile*& file, QLineEdit* edit)
{
    QString filename = edit->text();
    if (filename.isEmpty()) filename = deviceName.replace(QRegularExpression("[.,:=]"), "_");
    if (!filename.contains('.')) filename += ".dat";
    if (!filename.contains(QRegularExpression("[/\\\\]"))) filename = DefaultDirectory + filename;
    if (nullptr != file) delete file;
    file = new QFile(filename, this);
}

bool SoundRecordAndDrawControl::DetermineSampleTypeAndSize()
{
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    int deviceIndex = mInputSelectorBox->currentIndex();
	QAudioFormat format = deviceList[deviceIndex].preferredFormat();
	mSampleSize = format.bytesPerSample();
	mSampleType = format.sampleFormat();
    mSampleRate = format.sampleRate();
    mProcessedUSec = 0u;
    return true;
}

void SoundRecordAndDrawControl::StartRecording()
{
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    int deviceIndex = mInputSelectorBox->currentIndex();
    VerifyFileExists(mInputSelectorBox->currentText(), mOutputFile, mOutputFileNameEdit);
    if (!DetermineSampleTypeAndSize()) return;
    m_captureSession.audioInput()->setDevice(deviceList[deviceIndex]);
	m_audioRecorder->setOutputLocation(QUrl::fromLocalFile(mOutputFileNameEdit->text()));
	QMediaFormat format;
    format.setFileFormat(mFileFormatBox->currentData().value<QMediaFormat::FileFormat>());
    format.setAudioCodec(mCodecBox->currentData().value<QMediaFormat::AudioCodec>());
    m_audioRecorder->setMediaFormat(format);
    m_audioRecorder->setAudioSampleRate(mSampleRate);
    m_audioRecorder->setAudioChannelCount(1);
    m_audioRecorder->setEncodingMode(QMediaRecorder::ConstantBitRateEncoding);
	mOutputFile->open(QIODevice::WriteOnly);
    writeRST();
	m_audioRecorder->record();
    mInputSelectorBox->setEnabled(false);
    mStartButton->setEnabled(false);
    mStopButton->setEnabled(true);
}

void SoundRecordAndDrawControl::writeRST()
{
    QDataStream outputStream(mOutputFile);
    outputStream << RST.toUtf8() << static_cast<qint32>(mSampleRate) << static_cast<qint32>(mSampleSize) << static_cast<quint8>(mSampleType);
}

void SoundRecordAndDrawControl::Stop()
{
    mProcessedUSec = m_audioRecorder->duration() * 1000;
    m_audioRecorder->stop();
    mOutputFile->close();
    mInputSelectorBox->setEnabled(true);
    mStartButton->setEnabled(true);
    mStopButton->setEnabled(false);
}

void SoundRecordAndDrawControl::SplitFileIntoPackets()
{
    QString filename = mInputFileNameEdit->text();
    if (filename.isEmpty()) return;
    if (nullptr != mInputFile) delete mInputFile;
    mInputFile = new QFile(filename);
    if (!mInputFile->exists()) return;
    int allBytes = mInputFile->size(), fileIndex, currentOffset, packetSize = mPacketSizeEdit->text().toInt() * 1000000, n = filename.lastIndexOf('.'), headerSize = 16;
    if (0 == packetSize) return;
    mInputFile->open(QIODevice::ReadOnly);
    QDataStream inputStream(mInputFile);
    char header[headerSize];
    int bytesRead = inputStream.readRawData(header, headerSize);
    if (headerSize != bytesRead)
    {
        QMessageBox::information(this, "DrawSound", "The input file could not be read!");
        return;
    }
    if (RST == QString::fromUtf8(header + 4, 3)) mSampleSize = qFromBigEndian<qint32>(header + 11);
    else headerSize = 0;
    if (0 == mSampleSize && !DetermineSampleTypeAndSize()) return;
    packetSize = mSampleSize * (packetSize / mSampleSize);
    QString fileEnding = filename.right(filename.length() - n), basicFilename = filename.left(n);
    char* buffer = new char[packetSize];
    if (nullptr != mOutputFile)
    {
        delete mOutputFile;
        mOutputFile = nullptr;
    }
    for (fileIndex = 0, currentOffset = 0; currentOffset < allBytes; ++fileIndex)
    {
        QString splitFilename = basicFilename + QString::number(fileIndex) + fileEnding;
        mOutputFile = new QFile(splitFilename);
        mOutputFile->open(QIODevice::WriteOnly);
        int bytesRead = inputStream.readRawData(buffer, packetSize);
        if (-1 == bytesRead)
        {
            QMessageBox::information(this, "DrawSound", "The data could not be read completely!");
            break;
        }
        if (0 == fileIndex || 0 < headerSize)
        {
            currentOffset += bytesRead;
            int bytesWritten = mOutputFile->write(header, 12);
            if (bytesWritten < 12)
            {
                QMessageBox::information(this, "DrawSound", "Not all data could be written!");
                break;
            }
        }
        currentOffset += bytesRead;
        int bytesWritten = mOutputFile->write(buffer, bytesRead);
        if (bytesWritten < bytesRead)
        {
            QMessageBox::information(this, "DrawSound", "Not all data could be written!");
            break;
        }
        delete mOutputFile;
        mOutputFile = nullptr;
    }
    delete[] buffer;
}

void SoundRecordAndDrawControl::Save(const char *const inputData, const int nBytes)
{
    QString filename = mOutputFileNameEdit->text();
    if (filename.isEmpty())
    {
        filename = mInputFileNameEdit->text();
        mOutputFileNameEdit->setText(filename);
    }
    if (nullptr != mOutputFile) delete mOutputFile;
    mOutputFile = new QFile(filename);
    mOutputFile->open(QIODevice::WriteOnly);
    writeRST();
    int bytesWritten = mOutputFile->write(inputData, nBytes);
    if (bytesWritten < nBytes)
    {
        QMessageBox::information(this, "DrawSound", "Not all data could be written!");
    }
}

void SoundRecordAndDrawControl::createDecoder()
{
    if (nullptr == mDecoder)
    {
        mDecoder = new QAudioDecoder;
        connect(mDecoder, SIGNAL(finished()), this, SLOT(ReadyToDraw()));
        connect(mDecoder, SIGNAL(error(QAudioDecoder::Error)), this, SLOT(Error(QAudioDecoder::Error)));
        connect(mDecoder, SIGNAL(bufferReady()), this, SLOT(BufferReady()));
    }
}

void SoundRecordAndDrawControl::Decode()
{
    QString inputFileName = mInputFileNameEdit->text();
    if (inputFileName.isEmpty())
    {
        QMessageBox::information(this, "DrawSound", "Please provide an existing input file.");
        return;
    }
    if (nullptr != mInputFile) delete mInputFile;
    mInputFile = new QFile(inputFileName);
    if (!mInputFile->exists())
    {
        QMessageBox::information(this, "DrawSound", "The given input file does not exist!");
        return;
    }
    QString outputFileName = mOutputFileNameEdit->text();
    if (outputFileName.isEmpty())
    {
        int pointIndex = inputFileName.lastIndexOf('.');
        if (pointIndex < 1) outputFileName = inputFileName + ".dat";
        else outputFileName = inputFileName.left(pointIndex) + ".dat";
        mOutputFileNameEdit->setText(outputFileName);
    }
    if (nullptr != mOutputFile) delete mOutputFile;
    mOutputFile = new QFile(outputFileName);
    mInputFile->open(QIODevice::ReadOnly);
    mOutputFile->open(QIODevice::WriteOnly);
    createDecoder();
    mDecoder->setSourceDevice(mInputFile);
    mDecoder->start();
}

void SoundRecordAndDrawControl::Draw()
{
    VerifyFileExists(mInputSelectorBox->currentText(), mInputFile, mInputFileNameEdit);
    mInputFile->open(QIODevice::ReadOnly);
    if (mInputFile->fileName().right(4) == ".mp3" || mInputFile->fileName().right(4) == ".m4a")
    {
        createDecoder();
        mDecoder->setSourceDevice(mInputFile);
        mDecodeBuffer.clear();
        mProcessedUSec = 0;
        mDecodingFor = DF_Draw;
        mDecoder->start();
    }
    else
    {
        int nBytes = mInputFile->size();
        if (0 < nBytes)
        {
            QDataStream stream(mInputFile);
            char* inputData = new char[nBytes];
            nBytes = stream.readRawData(inputData, nBytes);
            int offset = 0;
            if (RST == QString::fromUtf8(inputData + 4, 3))
            {
                mSampleRate = qFromBigEndian<qint32>(inputData + 7);
                mSampleSize = qFromBigEndian<qint32>(inputData + 11);
                mSampleType = static_cast<QAudioFormat::SampleFormat>(inputData[15]);
                offset = 16;
            }
            mNumChannels = 1;
            if (0 == mSampleSize)
            {
                DefineSampleSizeDialog *dialog = new DefineSampleSizeDialog(this, inputData, nBytes);
                mMW->showMDIChild(dialog);
                return;
            }
            draw(inputData + offset, nBytes - offset);
            delete[] inputData;
        }
        else QMessageBox::information(this, "DrawSound", "With the selected input device no sound could be recorded!");
    }
}

void SoundRecordAndDrawControl::draw(const char* const inputData, const int nBytes)
{
    int nSamples = nBytes / (mSampleSize * mNumChannels) * 8, i, n, pos, p, sampleSizeBytes = mSampleSize / 8, frameSize = sampleSizeBytes * mNumChannels;
    if (nullptr != inputData && 0 < nBytes)
    {
        if (0u == mProcessedUSec) mProcessedUSec = static_cast<double>(nBytes) * 8000000 / (mSampleSize * mSampleRate * mNumChannels);
        double passedTime = 0.000001 * mProcessedUSec, currentTime;
        double timeStep = passedTime / nSamples;
        double ***data = Create(mNumChannels, nSamples, 2);
        for (i=pos=0, currentTime = 0.0; i < nSamples; ++i, currentTime += timeStep, pos += frameSize)
        {
            for (n=0; n < mNumChannels; ++n) data[n][i][0] = currentTime;
            switch (mSampleType)
            {
                case QAudioFormat::UnSignedInt:
                    switch (mSampleSize)
                    {
                        case 8:
                            for (n=0; n < mNumChannels; ++n) data[n][i][1] = static_cast<u_int8_t>(inputData[pos + n]);
                            break;
                        case 16:
                            for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes) data[n][i][1] = *reinterpret_cast<const u_int16_t*>(inputData + p);
                            break;
                        default:
                            for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes)
                            {
                                u_int32_t sample = 0;
                                memcpy(&sample, inputData + p, sampleSizeBytes);
                                data[n][i][1] = sample;
                                break;
                            }
                    }
                case QAudioFormat::Float:
                    if (mSampleSize < 32) for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes)
                    {
                        _Float32 sample = 0.0f;
                        memcpy(&sample, inputData + p, sampleSizeBytes);
                        data[n][i][1] = sample;
                    }
                    else if (mSampleSize == 32) for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes) data[n][i][1] = *reinterpret_cast<const _Float32*>(inputData + p);
                    else if (mSampleSize < 64) for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes)
                    {
                        double sample = 0.0;
                        memcpy(&sample, inputData + p, sampleSizeBytes);
                        data[n][i][1] = sample;
                    }
                    else if (mSampleSize == 64) for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes) data[n][i][1] = *reinterpret_cast<const double*>(inputData + p);
                    // printf("(%g|%g), ", data[i][0], data[i][1]);
                    break;
                case QAudioFormat::Unknown:
                case QAudioFormat::SignedInt:
                default:
                    switch (mSampleSize)
                    {
                        case 8:
                            for (n=0; n < mNumChannels; ++n) data[n][i][1] = inputData[i + n];
                            break;
                        case 16:
                            for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes) data[n][i][1] = *reinterpret_cast<const int16_t*>(inputData + p);
                            break;
                        default:
                            for (n=0, p=pos; n < mNumChannels; ++n, p += sampleSizeBytes)
                            {
                                int32_t sample = 0;
                                memcpy(&sample, inputData + p, sampleSizeBytes);
                                data[n][i][1] = sample;
                            }
                            break;
                    }
            }
        }
        SoundDrawWindow* window = new SoundWindow(this, mMW, mInputFileNameEdit->text(), mSampleRate);
        window->setData(data[0], nSamples);
        for (n=1; n < mNumChannels; ++n) window->addData(data[n], nSamples);
        delete data;
        mMW->showMDIChild(window);
    }
}

void SoundRecordAndDrawControl::Draw(const int sampleSize, const int sampleRate, const QAudioFormat::SampleType sampleType, const char *const inputData, const int nBytes)
{
    mSampleSize = sampleSize;
    mSampleRate = sampleRate;
    mSampleType = sampleType;
    draw(inputData, nBytes);
}

void SoundRecordAndDrawControl::ReadyToDraw()
{
    if (mDecodingFor == DF_Draw) draw(mDecodeBuffer.data(), mDecodeBuffer.size());
    else
    {
        delete mOutputFile;
        mOutputFile = nullptr;
    }
    mDecodingFor = DF_Nothing;
}

void SoundRecordAndDrawControl::Error(QAudioDecoder::Error)
{
    emit showMessage(DecodeError);
}

void SoundRecordAndDrawControl::BufferReady()
{
    QAudioBuffer buffer = mDecoder->read();
    if (0 == mProcessedUSec)
    {
        QAudioFormat format = buffer.format();
        mNumChannels = format.channelCount();
        if (mNumChannels > 1 && mDecodingFor != DF_Draw)
        {
            printf("Observed a channel count = %d. This is currently not supported!\n", format.channelCount());
            emit showMessage(ChannelCountLargerOne);
            mDecoder->stop();
            return;
        }
        mSampleRate = format.sampleRate();
        mSampleSize = format.sampleSize();
        mSampleType = format.sampleType();
        if (mDecodingFor == DF_File) writeRST();
    }
    if (mDecodingFor == DF_File)
    {
        static int callCounter = 0;
        printf("Number of Calls: %d\n", ++callCounter);
        int count = buffer.byteCount();
        int bytesWritten = mOutputFile->write(reinterpret_cast<const char*>(buffer.constData()), count);
        if (bytesWritten < count)
        {
            emit showMessage(NotAllDataCouldBeWritten);
            mDecoder->stop();
        }
    }
    else
    {
        mProcessedUSec += buffer.duration();
        mDecodeBuffer.append(reinterpret_cast<const char*>(buffer.constData()), buffer.byteCount());
    }
}

void SoundRecordAndDrawControl::ShowMessage(Message message)
{
    switch(message)
    {
        case ChannelCountLargerOne:
            QMessageBox::warning(this, "DrawSound", "Observed a channel count larger than one. This is currently not supported!");
            break;
        case NotAllDataCouldBeWritten:
            QMessageBox::warning(this, "DrawSound", "Error: Not all data could be written!");
            break;
        case DecodeError:
            QMessageBox::warning(this, "DrawSound", "Error: %s\n", mDecoder->errorString().toLatin1().data());
            break;
    }
}

void SoundRecordAndDrawControl::clearAssignmentData()
{
    if (nullptr != mAssignmentElements)
    {
        for (int n=0; n < mNumAssignmentElements; ++n) if (nullptr != mAssignmentElements[n].data) delete[] mAssignmentElements[n].data;
        delete[] mAssignmentElements;
    }
}

void SoundRecordAndDrawControl::InitializeAssignmentData(const int numElements, const int numOscillators)
{
    clearAssignmentData();
    mAssignmentElements = new AssignmentElement[numElements];
    for (int n=0; n < numElements; ++n)
    {
        mAssignmentElements[n].data = new double[numOscillators];
        for (int m=0; m < numOscillators; ++m) mAssignmentElements[n].data[m] = 0.0;
    }
}

SoundRecordAndDrawControl::AssignmentElement* SoundRecordAndDrawControl::GetAssignmentData(int& numElements) const
{
    numElements = mNumAssignmentElements;
    return mAssignmentElements;
}
