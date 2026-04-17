//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "definesamplesizedialog.h"
#include "recordanddrawControl.h"

#include <QComboBox>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
<<<<<<< HEAD
#include <QAudioDevice>
#include <QMediaDevices>


DefineSampleSizeDialog::DefineSampleSizeDialog(SoundRecordAndDrawControl* parent, const char *const inputData, const int nBytes) : QWidget(parent), mDeviceBox(new QComboBox(this)),  mSampleSizeBox(new QComboBox(this)),
    mSampleRateBox(new QComboBox(this)), mControl(parent), mInputData(inputData), mNBytes(nBytes)
{
    setAttribute(Qt::WA_DeleteOnClose);
    QGridLayout* layout = new QGridLayout(this);
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    for (QAudioDevice info : deviceList) mDeviceBox->addItem(info.description());
=======
#include <QAudioDeviceInfo>


DefineSampleSizeDialog::DefineSampleSizeDialog(SoundRecordAndDrawControl* parent, const char *const inputData, const int nBytes) : QWidget(parent), mDeviceBox(new QComboBox(this)), mSampleSizeBox(new QComboBox(this)),
    mSampleRateBox(new QComboBox(this)), mSampleTypeBox(new QComboBox(this)), mControl(parent), mInputData(inputData), mNBytes(nBytes)
{
    setAttribute(Qt::WA_DeleteOnClose);
    QGridLayout* layout = new QGridLayout(this);
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    for (QAudioDeviceInfo info : deviceList) mDeviceBox->addItem(info.deviceName());
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    mDeviceBox->setEditable(false);
    layout->addWidget(new QLabel("Device:", this), 0, 0);
    layout->addWidget(mDeviceBox, 0, 1);
    mSampleSizeBox->setEditable(false);
<<<<<<< HEAD
    layout->addWidget(new QLabel("Sample format:", this), 1, 0);
    layout->addWidget(mSampleSizeBox, 1, 1);
    layout->addWidget(new QLabel("Sample rate:", this), 2, 0);
    layout->addWidget(mSampleRateBox, 2, 1);
    layout->setRowMinimumHeight(3, 20);
    QPushButton* DrawB = new QPushButton("Draw", this), *SaveB = new QPushButton("Save", this), *CloseB = new QPushButton("Close", this);
    QGridLayout* BL = new QGridLayout;
    layout->addLayout(BL, 4, 0, 1, 2);
=======
    layout->addWidget(new QLabel("Sample size:", this), 1, 0);
    layout->addWidget(mSampleSizeBox, 1, 1);
    mSampleRateBox->setEditable(false);
    layout->addWidget(new QLabel("Sample rate:", this), 2, 0);
    layout->addWidget(mSampleRateBox, 2, 1);
    mSampleTypeBox->setEditable(false);
    layout->addWidget(new QLabel("Sample type:", this), 3, 0);
    layout->addWidget(mSampleTypeBox, 3, 1);
    layout->setRowMinimumHeight(4, 20);
    QPushButton* DrawB = new QPushButton("Draw", this), *SaveB = new QPushButton("Save", this), *CloseB = new QPushButton("Close", this);
    QGridLayout* BL = new QGridLayout;
    layout->addLayout(BL, 5, 0, 1, 2);
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    BL->addWidget(DrawB, 0, 0);
    BL->addWidget(SaveB, 0, 1);
    BL->addWidget(CloseB, 0, 2);
    connect(mDeviceBox, SIGNAL(currentIndexChanged(int)), this, SLOT(DeviceChanged(int)));
    connect(DrawB, SIGNAL(clicked()), this, SLOT(Draw()));
    connect(SaveB, SIGNAL(clicked()), this, SLOT(Save()));
    connect(CloseB, SIGNAL(clicked()), this, SLOT(close()));
    DeviceChanged(0);
}

DefineSampleSizeDialog::~DefineSampleSizeDialog() noexcept
{
    delete[] mInputData;
}

void DefineSampleSizeDialog::DeviceChanged(int index)
{
<<<<<<< HEAD
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    QList<QAudioFormat::SampleFormat> sizes = deviceList[index].supportedSampleFormats();
    mSampleSizeBox->clear();
	for (QAudioFormat::SampleFormat format : sizes)
	{
		switch (format)
		{
			case QAudioFormat::Unknown:
                mSampleSizeBox->addItem("unknown");
				break;
			case QAudioFormat::UInt8:
				mSampleSizeBox->addItem("UInt8");
				break;
			case QAudioFormat::Int16:
				mSampleSizeBox->addItem("Int16");
				break;
			case QAudioFormat::Int32:
				mSampleSizeBox->addItem("Int32");
				break;
			case QAudioFormat::Float:
				mSampleSizeBox->addItem("Float");
				break;
			default:
				// do nothing
				break;
		}
	}
	constexpr auto allSamplingRates = std::array{
        8000,  11025, 12000, 16000, 22050,  24000,  32000,  44100,
        48000, 64000, 88200, 96000, 128000, 176400, 192000,
    };
	mSampleRateBox->clear();
	int minSampleRate = deviceList[index].minimumSampleRate(), maxSampleRate = deviceList[index].maximumSampleRate();
	for (auto rate : allSamplingRates) if (rate >= minSampleRate && rate <= maxSampleRate) mSampleRateBox->addItem(QString::number(rate), rate);
=======
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    QList<int> sizes = deviceList[index].supportedSampleSizes();
    mSampleSizeBox->clear();
    for (int size : sizes) mSampleSizeBox->addItem(QString::number(size));
    mSampleRateBox->clear();
    QList<int> rates = deviceList[index].supportedSampleRates();
    for (int rate : rates) mSampleRateBox->addItem(QString::number(rate));
    mSampleTypeBox->clear();
    QList<QAudioFormat::SampleType> types = deviceList[index].supportedSampleTypes();
    for (auto type : types)
    {
        switch(type)
        {
            case QAudioFormat::Unknown:
                mSampleTypeBox->addItem("unknown");
                break;
            case QAudioFormat::SignedInt:
                mSampleTypeBox->addItem("signed int");
                break;
            case QAudioFormat::UnSignedInt:
                mSampleTypeBox->addItem("unsigned int");
                break;
            case QAudioFormat::Float:
                mSampleTypeBox->addItem("float");
                break;
        }
    }
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
}

void DefineSampleSizeDialog::Draw()
{
<<<<<<< HEAD
    QList<QAudioDevice> deviceList = QMediaDevices::audioInputs();
    int index = mDeviceBox->currentIndex();
    QList<QAudioFormat::SampleFormat> sizes = deviceList[index].supportedSampleFormats();
	QAudioFormat format = deviceList[index].preferredFormat();
	QAudioFormat::SampleFormat sampleFormat = sizes[mSampleSizeBox->currentIndex()];
	format.setSampleFormat(sampleFormat);
	format.setChannelCount(1);	
    char* data = new char[mNBytes];
    memcpy(data, mInputData, mNBytes);
    mControl->Draw(format.bytesPerSample(), mSampleRateBox->currentText().toInt(), sampleFormat, data, mNBytes);
=======
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    int index = mDeviceBox->currentIndex();
    QList<int> sizes = deviceList[index].supportedSampleSizes();
    QList<int> rates = deviceList[index].supportedSampleRates();
    QList<QAudioFormat::SampleType> types = deviceList[index].supportedSampleTypes();
    char* data = new char[mNBytes];
    memcpy(data, mInputData, mNBytes);
    mControl->Draw(sizes[mSampleSizeBox->currentIndex()], rates[mSampleRateBox->currentIndex()], types[mSampleTypeBox->currentIndex()], data, mNBytes);
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
}

void DefineSampleSizeDialog::Save()
{
    mControl->Save(mInputData, mNBytes);
}
