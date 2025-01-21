#include "soundwindow.h"
#include "frequencywindow.h"
#include "datensatz.h"
#include "utils.h"
#include "fit.h"
#include "windowselectdialog.h"
#include "soundvector.h"
#include "soundmatrix.h"
#include "roundbuffer.h"
#include "boxfilterdialog.h"
#include "recordanddrawControl.h"
#include "soundmainwindow.h"
#include "maxbuffer.h"

#include <QAudioOutput>
#include <QAction>
#include <QMenu>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QPainter>


SoundWindow::SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate) : SoundDrawWindow(control, sampleRate, 1), mOutputDeviceBox(new QComboBox(this)),
    mAudioOutput(nullptr), mAudioInputDevice(nullptr), mFilename(filename), mLabelOrderFilename(DATA_DIRECTORY "/Labels/Label.index"), mAddLabelAct(new QAction("Add label...", this)),
    mSaveLabelsAct(new QAction("Save labels (...)", this)), mDeleteAct(new QAction("Delete", this)), mPlayState(PSStopPlaying), mMinLabelWidth(-1.0)
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
    QAction *playAct = new QAction("Play", this), *FFTAct = new QAction("FFT", this), *FastLabelingAct = new QAction("Fast labeling", this), *LoadLabelsAct = new QAction("Load labels (...)", this);
    QAction *WriteAnnInputAct = new QAction("Write ANN input...", this), *ReadAnnOutputAct = new QAction("Read ANN output...", this), *ApplyBoxFilterAct = new QAction("Apply box filter...", this);
    QAction *ApplyDiffMaxTransfoAct = new QAction("Apply diff max transformation", this);
    FFTAct->setCheckable(true);
    mPopupMenu->addAction(playAct);
    mPopupMenu->addAction(FFTAct);
    mPopupMenu->addAction(ApplyBoxFilterAct);
    mPopupMenu->addAction(ApplyDiffMaxTransfoAct);
    FastLabelingAct->setCheckable(true);
    mPopupMenu->addAction(FastLabelingAct);
    mPopupMenu->addSeparator();
    mPopupMenu->addAction(mAddLabelAct);
    mPopupMenu->addAction(LoadLabelsAct);
    mPopupMenu->addAction(mSaveLabelsAct);
    mPopupMenu->addSeparator();
    mPopupMenu->addAction(WriteAnnInputAct);
    mPopupMenu->addAction(ReadAnnOutputAct);
    mPopupMenu->addSeparator();
    mPopupMenu->addAction(mDeleteAct);
    connect(playAct, SIGNAL(triggered()), this, SLOT(play()));
    connect(FFTAct, SIGNAL(toggled(bool)), this, SLOT(FFTActTriggered(bool)));
    connect(ApplyBoxFilterAct, SIGNAL(triggered()), this, SLOT(ApplyBoxFilter()));
    connect(ApplyDiffMaxTransfoAct, SIGNAL(triggered()), this, SLOT(ApplyDiffMaxTransfo()));
    connect(FastLabelingAct, SIGNAL(toggled(bool)), this, SLOT(setFastAssignmentMode(bool)));
    connect(mAddLabelAct, SIGNAL(triggered()), this, SLOT(AddLabel()));
    connect(LoadLabelsAct, SIGNAL(triggered()), this, SLOT(LoadLabels()));
    connect(mSaveLabelsAct, SIGNAL(triggered()), this, SLOT(SaveLabels()));
    connect(WriteAnnInputAct, SIGNAL(triggered()), this, SLOT(WriteAnnInput()));
    connect(mDeleteAct, SIGNAL(triggered()), this, SLOT(Delete()));
    connect(ReadAnnOutputAct, SIGNAL(triggered()), this, SLOT(ReadAndVerifyAnnOutput()));
    QFile labelOrderFile(mLabelOrderFilename);
    labelOrderFile.open(QIODevice::ReadOnly);
    mLabelOrder = QString(labelOrderFile.readAll()).split('\t');
    if (!mLabelOrder.empty())
    {
        mMinLabelWidth = mLabelOrder[0].toDouble();
        mLabelOrder.pop_front();
    }
}

SoundWindow::~SoundWindow()
{
    if (nullptr != mAudioOutput) delete mAudioOutput;
    else if (nullptr != mAudioInputDevice) delete mAudioInputDevice;
}

void SoundWindow::play()
{
    mPlayState = PSPlayOnce;
    startPlaying();
}

void SoundWindow::setFastAssignmentMode(bool enable)
{
    if (enable)
    {
        mMode = MFastLabeling;
        mPlayState = PSPlayContinuously;
        int selWidth = mMinLabelWidth * XSF;
        if (nullptr == mSelectionRect)
        {
            QPoint mousePosition(mapFromGlobal(QCursor::pos()));
            int x = mousePosition.x() - 0.5 * selWidth;
            if (x < ScaleYWidth) x = ScaleYWidth;
            else if (x + selWidth > width()) x = width() - selWidth;
            QRect newSelection(x, 0, selWidth, height() - ScaleXHeight);
            SelectionChanged(&newSelection);
        }
        else mSelectionRect->setWidth(selWidth);
        startPlaying();
    }
    else
    {
        mMode = MNormal;
        mPlayState = PSStopPlaying;
    }
}

void SoundWindow::startPlaying()
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
    else if (nullptr != mAudioInputDevice) delete mAudioInputDevice;
    mAudioOutput = new QAudioOutput(deviceList[mOutputDeviceBox->currentIndex()], format, this);
    mAudioInputDevice = mAudioOutput->start();
    if (mPlayState == PSPlayContinuously) connect(mAudioOutput, SIGNAL(notify()), this, SLOT(continuePlaying()));
    continuePlaying();
}

void SoundWindow::continuePlaying()
{
    switch(mPlayState)
    {
        case PSPlayOnce:
            mPlayState = PSStopPlaying;
            break;
        case PSPlayContinuously:
            break;
        case PSStopPlaying:
            return;
    }
    float* data;
    int length = getSoundData(&data);
    mAudioInputDevice->write(reinterpret_cast<char*>(data), 4*length);
    if (mPlayState == PSPlayContinuously) mAudioOutput->setNotifyInterval(static_cast<int>((static_cast<double>(length) / mSampleRate) * 1000));
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

int SoundWindow::getSoundData(double ** data, const int labelIndex, FFTSelection fftSelection) const
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
        double **realFFTData, **imaginaryFFTData;
        int FFTLength;
        getFFTData(FSForFTT, -1, FFTLength, realFFTData, imaginaryFFTData);
        if (nullptr == mFFTWindow)
        {
            mFFTWindow = new FrequencyWindow(mControl, mSampleRate);
            mControl->GetMW()->showMDIChild(mFFTWindow);
        }
        else mFFTWindow->clear();
        mFFTWindow->setData(realFFTData, FFTLength);
        mFFTWindow->addData(imaginaryFFTData, FFTLength);
        if (!mFFTWindow->isVisible()) mFFTWindow->show();
        Destroy(realFFTData, FFTLength);
        Destroy(imaginaryFFTData, FFTLength);
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
            double **realFFTData, **imaginaryFFTData;
            int FFTLength;
            getFFTData(FSForSelectedFTT, n, FFTLength, realFFTData, imaginaryFFTData);
            DiagWindow* window = new DiagWindow;
            window->setWindowTitle(mLabels[n].phoneme);
            window->setUnits("frequency [Hz]", "intensity");
            window->setData(realFFTData, FFTLength);
            window->addData(imaginaryFFTData, FFTLength);
            mControl->GetMW()->showMDIChild(window);
            Destroy(realFFTData, FFTLength);
            Destroy(imaginaryFFTData, FFTLength);
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
    addLabel(dialog.GetName());
}

void SoundWindow::addLabel(const QString name)
{
    Label newLabel;
    newLabel.phoneme = name;
    newLabel.rect = *mSelectionRect;
    if (mSelectionRect->width() < mMinLabelWidth) mMinLabelWidth = mSelectionRect->width();
    newLabel.index = estimateLabelIndex(newLabel.phoneme);
    mLabels.push_back(newLabel);
    if (mMode != MFastLabeling)
    {
        delete mSelectionRect;
        mSelectionRect = nullptr;
    }
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
        stream << label.phoneme << '[' << label.index << ']' << '(' << QString::number(label.rect.left()).replace(',', '.') << ", " << QString::number(label.rect.top()).replace(',', '.') << ", "
               << QString::number(label.rect.right()).replace(',', '.') << ", " << QString::number(label.rect.bottom()).replace(',', '.') << ")\n";
    QFile labelOrderFile(mLabelOrderFilename);
    labelOrderFile.open(QIODevice::WriteOnly);
    labelOrderFile.write((QString::number(mMinLabelWidth, 'f', 15) + '\t' + mLabelOrder.join('\t')).toLatin1());
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
        size_t indexIndexLeft(line.indexOf('[')), indexIndexRight(line.indexOf(']'));
        size_t indexLeft(line.indexOf('(')), indexRight(line.indexOf(')'));
        if (0 < indexLeft && indexLeft < indexRight)
        {
            QStringList list = line.mid(indexLeft + 1, indexRight - indexLeft - 1).split(',');
            if (list.size() == 4)
            {
                Label label;
                if (indexIndexLeft > 0 && indexIndexLeft < indexIndexRight - 1) label.index = line.mid(indexIndexLeft + 1, indexIndexRight - indexIndexLeft - 1).toInt();
                else indexIndexLeft = indexLeft;
                label.phoneme = line.left(indexIndexLeft).trimmed();
                label.rect.setCoords(list[0].toDouble(), list[1].toDouble(), list[2].toDouble(), list[3].toDouble());
                if (mMinLabelWidth < 0.0 || label.rect.width() < mMinLabelWidth) mMinLabelWidth = label.rect.width();
                label.index = estimateLabelIndex(label.phoneme);
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

void SoundWindow::mouseLeftClicked(QPoint* Position)
{
    SoundDrawWindow::mouseLeftClicked(Position);
    if (!mIsFFT || 0 == mLabels.size()) return;
    calcMinLabelWidth();
    showFFT();
}

void SoundWindow::CreateAnnInput(double ** data, int &FFTLength)
{
    calcMinLabelWidth();
    double **realFFTData, **imaginaryFFTData;
    int i;
    for (int n=0; n < mLabels.size(); ++n)
    {
        getFFTData(FSForSelectedFTT, n, FFTLength, realFFTData, imaginaryFFTData);
        data[n] = new double[FFTLength];
        IntensityMax IMax[10];
        for(int m=0; m < FFTLength; ++m)
        {
            data[n][m] = (realFFTData[m][1] * realFFTData[m][1] + imaginaryFFTData[m][1] * imaginaryFFTData[m][1]);
            if (m>1 && data[n][m-1] > data[n][m-2] && data[n][m-1] > data[n][m])
            {
                for (i = 8; i >= 0 && data[n][m-1] > IMax[i].I; --i) IMax[i+1] = IMax[i];
                if (data[n][m-1] > IMax[i+1].I)
                {
                    IMax[i+1].I = data[n][m-1];
                    IMax[i+1].F = m-1;
                }
            }
        }
        Destroy(realFFTData, FFTLength);
        Destroy(imaginaryFFTData, FFTLength);
        printf("%s ", mLabels[n].phoneme.toLatin1().data());
        for (i=9; i>=0; --i) printf("%d ", IMax[i].F);
        printf("\n");
    }
}

void SoundWindow::WriteAnnInput()
{
    if (0 == mLabels.size())
    {
        QMessageBox::information(this, "DrawSound", "No ANN input can be created because no labels are loaded.");
        return;
    }
    QString filename = QFileDialog::getSaveFileName(this, "Select filename", predictLabelFilename());
    if (filename.isEmpty()) return;
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QDataStream stream(&file);
    int FFTLength;
    double **data = new double*[mLabels.size()];
    CreateAnnInput(data, FFTLength);
    for (int n=0; n < mLabels.size(); ++n)
    {
        if (n==0) stream << static_cast<quint32>(mLabels.size()) << static_cast<quint32>(FFTLength);
        stream << static_cast<quint8>(mLabels[n].index);
        for(int m=0; m < FFTLength; ++m) stream << data[n][m];
    }
    Destroy(data, mLabels.size());
}

void SoundWindow::ReadAndVerifyAnnOutput()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select filename", predictLabelFilename());
    if (filename.isEmpty()) return;
    QFile file(filename);
    file.open(QIODevice::ReadOnly);
    QDataStream stream(&file);
    quint16 T1NR, T1NC, T2NR, T2NC;
    stream >> T1NR >> T1NC;
    SoundMatrix T1(T1NC, T1NR);
    for (int c=0; c < T1NC; ++c) for (int r=0; r < T1NR; ++r) stream >> T1[c][r];
    stream >> T2NR >> T2NC;
    SoundMatrix T2(T2NC, T2NR);
    for (int c=0; c < T2NC; ++c) for (int r=0; r < T2NR; ++r) stream >> T2[c][r];

    int FFTLength, correctPred=0;
    double **data = new double*[mLabels.size()];
    CreateAnnInput(data, FFTLength);
    for (int n=0; n < mLabels.size(); ++n)
    {
        SoundVector in(FFTLength + 1);
        in[0] = 1.0;
        for (int m=0; m < FFTLength; ++m) in[m+1] = data[n][m];
        SoundVector h1((T1*in).sigmoid());
        if (mLabels[n].index == (T2*h1).getIndexOfMax()) ++correctPred;
    }
    Destroy(data, mLabels.size());

    QMessageBox::information(this, "DrawSound", QString("%1 of %2 predicted correctly!").arg(correctPred).arg(mLabels.size()));
}

void SoundWindow::calcMinLabelWidth()
{
    double minSize = mLabels[0].rect.width();
    for (Label label : mLabels) if (label.rect.width() < minSize) minSize = label.rect.width();
    mSelectedFFtSize = getFFTWidth(minSize);
}

void SoundWindow::getFFTData(const SoundDrawWindow::FFTSelection selection, const int labelIndex, int& FFTLength, double **& realData, double **& imaginaryData) const
{
    double* data;
    int length = getSoundData(&data, labelIndex, selection);
    FFTLength = length + 1;
    realData = Create(FFTLength, 2);
    imaginaryData = Create(FFTLength, 2);
    calcFFT(data, length, 1.0 / mSampleRate, realData, imaginaryData);
    delete[] data;
}

void SoundWindow::keyPressed(QKeyEvent* K)
{
    if (mMode == MFastLabeling && K->text().length() == 1)
    {
        K->accept();
        mKeyText += K->text();
        QStringList results = mLabelOrder.filter(QRegularExpression("^" + mKeyText, QRegularExpression::CaseInsensitiveOption));
        if (results.size() == 1) addLabel(results[0]);
        else if (results.empty()) addLabel(mKeyText);
        else return;
        mKeyText.clear();
    }
}

void SoundWindow::setData(double ** Data, int numRows)
{
    analyzeData(Data, numRows);
    DiagWindow::setData(Data, numRows);
}

void SoundWindow::analyzeData(double **const Data, const int numRows)
{

}

void SoundWindow::ApplyBoxFilter()
{
    BoxFilterDialog dialog;
    if (dialog.exec() == QDialog::Rejected) return;
    const int filterRadius = dialog.getResult(), nData = Daten->GetDSL();
    double **filteredData = Create(nData, 2), Sum = 0.0;
    int l, c, r;
    for (r=0, c = -filterRadius, l = - 2 * filterRadius - 1; c < nData; ++r, ++c, ++l)
    {
        if (r < nData) Sum += Daten->GetValue(r, 1);
        if (c >= 0)
        {
            if (l >= 0) Sum -= Daten->GetValue(l, 1);
            filteredData[c][0] = Daten->GetValue(c, 0);
            filteredData[c][1] = Sum;
        }
    }

    SoundWindow* newWindow = new SoundWindow(mControl, mFilename, mSampleRate);
    newWindow->setWindowTitle(newWindow->windowTitle() + " BoxFilter " + QString::number(filterRadius));
    newWindow->setData(filteredData, nData);
    newWindow->mLabels = mLabels;
    mControl->GetMW()->showMDIChild(newWindow);
}

void SoundWindow::ApplyDiffMaxTransfo()
{
    const double radius_s = 0.005;
    const int radius = static_cast<int>(radius_s * mSampleRate), nData = Daten->GetDSL();
    double** transformedData = Create(nData, 2), t, a, step = 1.0 / mSampleRate;
    MaxBuffer maxbuffer(2 * radius_s);
    int c, r;
    for (r=0, c = -radius; c < nData; r++, c++)
    {
        if (r < nData)
        {
            t = Daten->GetValue(r, 0);
            a = abs(Daten->GetValue(r, 1));
        }
        else
        {
            t += step;
            a = 0.0;
        }
        if (c >= 0)
        {
            transformedData[c][0] = Daten->GetValue(c, 0);
            transformedData[c][1] = maxbuffer.newValue(t, a);
        }
        else maxbuffer.newValue(t, a);
    }
    SoundWindow* newWindow = new SoundWindow(mControl, mFilename, mSampleRate);
    newWindow->setWindowTitle(newWindow->windowTitle() + " DiffMaxTransformed");
    newWindow->setData(transformedData, nData);
    newWindow->mLabels = mLabels;
    mControl->GetMW()->showMDIChild(newWindow);
}
