//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
#include <QAudioDeviceInfo>


DefineSampleSizeDialog::DefineSampleSizeDialog(SoundRecordAndDrawControl* parent, const char *const inputData, const int nBytes) : QWidget(parent), mDeviceBox(new QComboBox(this)), mSampleSizeBox(new QComboBox(this)),
    mSampleRateBox(new QComboBox(this)), mSampleTypeBox(new QComboBox(this)), mControl(parent), mInputData(inputData), mNBytes(nBytes)
{
    setAttribute(Qt::WA_DeleteOnClose);
    QGridLayout* layout = new QGridLayout(this);
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    for (QAudioDeviceInfo info : deviceList) mDeviceBox->addItem(info.deviceName());
    mDeviceBox->setEditable(false);
    layout->addWidget(new QLabel("Device:", this), 0, 0);
    layout->addWidget(mDeviceBox, 0, 1);
    mSampleSizeBox->setEditable(false);
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
}

void DefineSampleSizeDialog::Draw()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioInput);
    int index = mDeviceBox->currentIndex();
    QList<int> sizes = deviceList[index].supportedSampleSizes();
    QList<int> rates = deviceList[index].supportedSampleRates();
    QList<QAudioFormat::SampleType> types = deviceList[index].supportedSampleTypes();
    char* data = new char[mNBytes];
    memcpy(data, mInputData, mNBytes);
    mControl->Draw(sizes[mSampleSizeBox->currentIndex()], rates[mSampleRateBox->currentIndex()], types[mSampleTypeBox->currentIndex()], data, mNBytes);
}

void DefineSampleSizeDialog::Save()
{
    mControl->Save(mInputData, mNBytes);
}
