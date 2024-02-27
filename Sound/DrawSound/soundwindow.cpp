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
    mAudioOutput(nullptr), mFilename(filename), mAddLabelAct(new QAction("Add label...", this)), mSaveLabelsAct(new QAction("Save labels (...)", this)), mDeleteAct(new QAction("Delete", this))
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
    mPopupMenu->addAction(mDeleteAct);
    connect(playAct, SIGNAL(triggered()), this, SLOT(Play()));
    connect(FFTAct, SIGNAL(toggled(bool)), this, SLOT(FFTActTriggered(bool)));
    connect(mAddLabelAct, SIGNAL(triggered()), this, SLOT(AddLabel()));
    connect(LoadLabelsAct, SIGNAL(triggered()), this, SLOT(LoadLabels()));
    connect(mSaveLabelsAct, SIGNAL(triggered()), this, SLOT(SaveLabels()));
    connect(mDeleteAct, SIGNAL(triggered()), this, SLOT(Delete()));
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
    int *xStart, *xStop, n, i=0, rLength = 0, nSelected = 0, j;
    if (nullptr != mSelectionRect)
    {
        nSelected = 1;
        xStart = new int[nSelected];
        xStop = new int[nSelected];
        rLength = getSoundDataRange(xStart[0], xStop[0]);
    }
    else
    {
        for (Label label : mLabels) if (label.isSelected) ++nSelected;
        if (0 < nSelected)
        {
            xStart = new int[nSelected];
            xStop = new int[nSelected];
            for (n=0; n < mLabels.size(); ++n) if (mLabels[n].isSelected)
            {
                rLength += getSoundDataRange(xStart[i], xStop[i], n, FSNotForFTT);
                ++i;
            }
        }
        else
        {
            nSelected = 1;
            xStart = new int[nSelected];
            xStop = new int[nSelected];
            xStart[0] = 0;
            rLength = Daten->GetDSL();
            xStop[0] = rLength - 1;
        }
    }
    *data = new float[rLength];
    for (i=j=0; j < nSelected; ++j) for (n = xStart[j]; n <= xStop[j]; ++n, ++i) (*data)[i] = static_cast<float>(Daten->GetValue(n, 1));
    delete[] xStart;
    delete[] xStop;
    return rLength;
}

int SoundWindow::getSoundData(double ** data, const int labelIndex, FFTSelection fftSelection)
{
    int xStart, xStop, n, i, rLength = getSoundDataRange(xStart, xStop, labelIndex, fftSelection);
    *data = new double[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = Daten->GetValue(n, 1);
    return rLength;
}

void SoundWindow::showFFT()
{
    if (nullptr != mSelectionRect)
    {
        double* data;
        int length = getSoundData(&data, -1, FSForFTT), FFTLength = length + 1;
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
    else if (mLabels.size() > 0)
    {
        int nSelected = 0;
        bool wasSelected = true;
        for (Label label : mLabels) if (label.isSelected) ++nSelected;
        if (0 == nSelected)
        {
            wasSelected = false;
            for (auto it = mLabels.begin(); it != mLabels.end(); ++it) it->isSelected = true;
            nSelected = mLabels.size();
        }
        for (int n=0; n < mLabels.size(); ++n) if (mLabels[n].isSelected)
        {
            double* data;
            int length = getSoundData(&data, n, FSForFTT), FFTLength = length + 1;
            double **realFFTData = Create(FFTLength, 2), **imaginaryFFTData = Create(FFTLength, 2);
            calcFFT(data, length, 1.0 / mSampleRate, realFFTData, imaginaryFFTData);
            DiagWindow* window = new DiagWindow;
            window->setWindowTitle(mLabels[n].phoneme);
            window->setUnits("frequency [Hz]", "intensity");
            window->setData(realFFTData, FFTLength);
            window->addData(imaginaryFFTData, FFTLength);
            window->show();
            Destroy(realFFTData, FFTLength);
            Destroy(imaginaryFFTData, FFTLength);
            delete[] data;
        }
        if (!wasSelected) for (auto it = mLabels.begin(); it != mLabels.end(); ++it) it->isSelected = false;
    }
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
    mLabelFont = P.font();
    mLabelFont.setPixelSize(24);
    for (Label label : mLabels)
    {
        P.setPen(label.isSelected ? QColor(255, 255, 0) : QColor(0, 200, 0));
        int bottom = YO - label.rect.bottom() * YSF, height = label.rect.height() * YSF;
        int left = label.rect.left() * XSF + XO, width = label.rect.width() * XSF;
        P.drawRect(left, bottom, width, height);
        QRect textRect = getLabelTextRect(label);
        WriteText(P, textRect.left(),  textRect.bottom(), label.phoneme, mLabelFont, 0);
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
    auto state = getBestMouseState(point);
    if (state.first != -1 && state.second != MSOutside && !mLabels[state.first].isSelected)
    {
        for (int i=0; i < mLabels.size(); ++i) mLabels[i].isSelected = (i == state.first);
        if (nullptr != mSelectionRect)
        {
            delete mSelectionRect;
            mSelectionRect = nullptr;
        }
    }
    mAddLabelAct->setEnabled(nullptr != mSelectionRect);
    mSaveLabelsAct->setEnabled(!isSaved() && 0 < mLabels.size());
    SoundDrawWindow::ShowPopupMenu(point);
}

void SoundWindow::Delete()
{
    if (nullptr != mSelectionRect) delete mSelectionRect;
    mSelectionRect = nullptr;
    std::vector<Label> tempLabels;
    for (Label label : mLabels) if (!label.isSelected) tempLabels.push_back(label);
    mLabels = tempLabels;
    Paint();
}
