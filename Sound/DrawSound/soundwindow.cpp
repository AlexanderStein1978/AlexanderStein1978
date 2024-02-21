#include "soundwindow.h"
#include "frequencywindow.h"
#include "datensatz.h"
#include "utils.h"
#include "fit.h"
#include "windowselectdialog.h"

#include <QAudioOutput>
#include <QAction>
#include <QMenu>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QPainter>


SoundWindow::SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate) : SoundDrawWindow(control, sampleRate, 1), mOutputDeviceBox(new QComboBox(this)),
    mAudioOutput(nullptr), mFilename(filename), mAddLabelAct(new QAction("Add label...", this)), mSaveLabelsAct(new QAction("Save labels (...)", this))
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    for (QAudioDeviceInfo info : deviceList) mOutputDeviceBox->addItem(info.deviceName());
    mOutputDeviceBox->setEditable(false);
    QGridLayout *Layout = new QGridLayout;
    Layout->addWidget(new QLabel("Sound output device:", this), 0, 0);
    Layout->addWidget(mOutputDeviceBox, 0, 1);
    SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
    setUnits("time [s]", "intensity");
    setWindowTitle("Draw sound: " + filename);
    QAction *playAct = new QAction("Play", this), *FFTAct = new QAction("FFT", this), *LoadLabelsAct = new QAction("Load labels (...)", this);
    FFTAct->setCheckable(true);
    mPopupMenu->addAction(playAct);
    mPopupMenu->addAction(FFTAct);
    mPopupMenu->addSeparator();
    mPopupMenu->addAction(mAddLabelAct);
    mPopupMenu->addAction(LoadLabelsAct);
    mPopupMenu->addAction(mSaveLabelsAct);
    mPopupMenu->addSeparator();
    connect(playAct, SIGNAL(triggered()), this, SLOT(Play()));
    connect(FFTAct, SIGNAL(toggled(bool)), this, SLOT(FFTActTriggered(bool)));
    connect(mAddLabelAct, SIGNAL(triggered()), this, SLOT(AddLabel()));
    connect(LoadLabelsAct, SIGNAL(triggered()), this, SLOT(LoadLabels()));
    connect(mSaveLabelsAct, SIGNAL(triggered()), this, SLOT(SaveLabels()));
}


SoundWindow::~SoundWindow()
{
    if (nullptr != mAudioOutput) delete mAudioOutput;
}

void SoundWindow::Play()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    QAudioFormat format;
    format.setSampleRate(mSampleRate);
    format.setChannelCount(1);
    format.setSampleSize(32);
    format.setCodec("audio/pcm");
    format.setByteOrder(QAudioFormat::LittleEndian);
    format.setSampleType(QAudioFormat::Float);
    if (nullptr != mAudioOutput) delete mAudioOutput;
    mAudioOutput = new QAudioOutput(deviceList[mOutputDeviceBox->currentIndex()], format, this);
    QIODevice* inputDevice = mAudioOutput->start();
    float* data;
    int length = getSoundData(&data);
    inputDevice->write(reinterpret_cast<char*>(data), 4*length);
    delete[] data;
}

int SoundWindow::getSoundData(float ** data)
{
    int xStart, xStop, n, i, rLength = getSoundDataRange(xStart, xStop);
    *data = new float[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = static_cast<float>(Daten->GetValue(n, 1));
    return rLength;
}

int SoundWindow::getSoundData(double ** data)
{
    int xStart, xStop, n, i, rLength = getSoundDataRange(xStart, xStop);
    *data = new double[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = Daten->GetValue(n, 1);
    return rLength;
}

void SoundWindow::showFFT()
{
    double* data;
    int length = getSoundData(&data), FFTLength = length + 1;
    double **realFFTData = Create(FFTLength, 2), **imaginaryFFTData = Create(FFTLength, 2);
    calcFFT(data, length, 1.0 / mSampleRate, realFFTData, imaginaryFFTData);
    if (nullptr == mFFTWindow) mFFTWindow = new FrequencyWindow(mControl, mSampleRate);
    else mFFTWindow->clear();
    mFFTWindow->setData(realFFTData, FFTLength);
    mFFTWindow->addData(imaginaryFFTData, FFTLength);
    if (!mFFTWindow->isVisible()) mFFTWindow->show();
    Destroy(realFFTData, FFTLength);
    Destroy(imaginaryFFTData, FFTLength);
    delete[] data;
}

void SoundWindow::FFTActTriggered(bool checked)
{
    mIsFFT = checked;
    if (checked && nullptr != mSelectionRect)
    {
        mSelectionRect->setWidth(getFFTWidth(mSelectionRect->width()));
        showFFT();
        Paint();
    }
}

void SoundWindow::WriteToFile()
{
    QString filename = QFileDialog::getSaveFileName(this, "Select filename to save", DATA_DIRECTORY  "/Chars");
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    float* data;
    int length = getSoundData(&data);
    file.write(reinterpret_cast<char*>(data), 4*length);
    delete[] data;
}

void SoundWindow::AddLabel()
{
    NameSelectionDialog dialog;
    if (dialog.exec() == QDialog::Rejected) return;
    Label newLabel;
    newLabel.phoneme = dialog.GetName();
    newLabel.rect = *mSelectionRect;
    mLabels.push_back(newLabel);
    delete mSelectionRect;
    mSelectionRect = nullptr;
    ensureMouseShape(Qt::ArrowCursor);
    mMouseState = mMoveState = MSOutside;
    Changed();
    Paint();
}

QString SoundWindow::predictLabelFilename()
{
    size_t index;
    return QString(DATA_DIRECTORY "/Labels/") + ((index = mFilename.indexOf('.')) > 0 ? mFilename.left(index) : mFilename) + ".label";
}

void SoundWindow::SaveLabels()
{
    if (mLabelFilename.isEmpty())
    {
        mLabelFilename = QFileDialog::getSaveFileName(this, "Select filename", predictLabelFilename());
        if (mLabelFilename.isEmpty()) return;
    }
    QFile file(mLabelFilename);
    if (!file.open(QIODevice::WriteOnly))
    {
        QMessageBox::information(this, "DrawSound", "An error happened during opening file for writing!");
        return;
    }
    QTextStream stream(&file);
    stream << "Sound file: " << mFilename << '\n' << "Phoneme(left, top, right, bottom)\n";
    for (Label label : mLabels)
        stream << label.phoneme << '(' << QString::number(label.rect.left()).replace(',', '.') << ", " << QString::number(label.rect.top()).replace(',', '.') << ", "
               << QString::number(label.rect.right()).replace(',', '.') << ", " << QString::number(label.rect.bottom()).replace(',', '.') << ")\n";
    Saved();
}

void SoundWindow::LoadLabels()
{
    QString filename(mLabelFilename.isEmpty() || !QFile::exists(mLabelFilename) ? predictLabelFilename() : mLabelFilename);
    if (!QFile::exists(filename))
    {
        filename = QFileDialog::getOpenFileName(this, "Select filename", filename);
        if (filename.isEmpty()) return;
    }
    mLabels.clear();
    QFile file(filename);
    file.open(QIODevice::ReadOnly);
    file.readLine();
    file.readLine();
    while(!file.atEnd())
    {
        QString line(file.readLine());
        size_t indexLeft(line.indexOf('(')), indexRight(line.indexOf(')'));
        if (0 < indexLeft && indexLeft < indexRight)
        {
            QStringList list = line.mid(indexLeft + 1, indexRight - indexLeft - 1).split(',');
            if (list.size() == 4)
            {
                Label label;
                label.phoneme = line.left(indexLeft).trimmed();
                label.rect.setCoords(list[0].toDouble(), list[1].toDouble(), list[2].toDouble(), list[3].toDouble());
                mLabels.push_back(label);
            }
        }
    }
    Paint();
    Saved();
}

void SoundWindow::PSpektrum(QPainter& P, const QRect& A, bool PrintFN)
{
    SoundDrawWindow::PSpektrum(P, A, PrintFN);
    P.setPen(QColor(0, 200, 0));
    QFont font = P.font();
    font.setPixelSize(24);
    for (Label label : mLabels)
    {
        int bottom = YO - label.rect.bottom() * YSF, height = label.rect.height() * YSF;
        int left = label.rect.left() * XSF + XO, width = label.rect.width() * XSF;
        P.drawRect(left, bottom, width, height);
        WriteText(P, left + (abs(width) - TextWidth(font, label.phoneme)) / 2,  bottom - abs(height) - TextHeight(font, label.phoneme) / 4, label.phoneme, font, 0);
    }
}

void SoundWindow::closeEvent(QCloseEvent* i_event)
{
    if (isSaved()) i_event->accept();
    else
    {
        QMessageBox::StandardButton button = QMessageBox::question(this, "Labels not saved!", "Save labels now?", QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        if (button == QMessageBox::Yes)
        {
            SaveLabels();
            if (!isSaved()) closeEvent(i_event);
            else i_event->accept();
        }
        else if (button == QMessageBox::No) i_event->accept();
        else i_event->ignore();
    }
}

void SoundWindow::ShowPopupMenu(const QPoint& point)
{
    mAddLabelAct->setEnabled(nullptr != mSelectionRect);
    mSaveLabelsAct->setEnabled(!isSaved() && 0 < mLabels.size());
    SoundDrawWindow::ShowPopupMenu(point);
}
