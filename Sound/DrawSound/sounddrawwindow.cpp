#include "sounddrawwindow.h"

#include "datensatz.h"

#include <QPainter>
#include <QToolButton>
#include <QComboBox>
#include <QAudioOutput>
#include <QMenu>
#include <QFileDialog>


SoundDrawWindow::SoundDrawWindow(const QString& filename, const int sampleRate) : DiagWindow(SimpleDiagWindow, nullptr, "Data files (*.dat)", ".dat", 1),
    mPopupMenu(new QMenu(this)), mOutputDeviceBox(new QComboBox(this)), mAudioOutput(nullptr), mSampleRate(sampleRate)
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    for (QAudioDeviceInfo info : deviceList) mOutputDeviceBox->addItem(info.deviceName());
    QGridLayout *Layout = new QGridLayout;
    Layout->addWidget(new QLabel("Sound output device:", this), 0, 0);
    Layout->addWidget(mOutputDeviceBox, 0, 1);
    SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
    setUnits("time [s]", "intensity");
    setWindowTitle("Draw sound: " + filename);
	QAction *playAct = new QAction("Play", this), *writeAct = new QAction("Write to file...", this);
    mPopupMenu->addAction(playAct);
    mPopupMenu->addAction(writeAct);
	connect(playAct, SIGNAL(triggered()), this, SLOT(Play()));
    connect(writeAct, SIGNAL(triggered()), this, SLOT(WriteToFile()));
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(SelectionChanged(QRect*)));
    connect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
    connect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
    connect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
}

SoundDrawWindow::~SoundDrawWindow() noexcept
{
    if (nullptr != mSelectionRect) delete mSelectionRect;
    if (nullptr != mAudioOutput) delete mAudioOutput;
}

void SoundDrawWindow::SelectionChanged(QRect* MarkedArea)
{
    if (!ZoomB->isChecked())
    {
        if (nullptr == mSelectionRect) mSelectionRect = new QRectF;
        mSelectionRect->setCoords(double(MarkedArea->left() - XO) / XSF, -double(MarkedArea->top() - YO) / YSF, double(MarkedArea->right() - XO) / XSF, -double(MarkedArea->bottom() - YO) / YSF);
        Paint();
    }
}

void SoundDrawWindow::PSpektrum(QPainter& P, const QRect& A, bool PrintFN)
{
    DiagWindow::PSpektrum(P, A, PrintFN);
    if (nullptr != mSelectionRect)
        P.fillRect(mSelectionRect->left() * XSF + XO, YO - mSelectionRect->bottom() * YSF, mSelectionRect->width() * XSF, mSelectionRect->height() * YSF, QColor(127, 127, 127, 127));
}

void SoundDrawWindow::ensureMouseShape(const Qt::CursorShape shape)
{
    QCursor mousePointer = Bild->cursor();
    Qt::CursorShape currentShape = mousePointer.shape();
    if (currentShape != shape)
    {
        mousePointer.setShape(shape);
        Bild->setCursor(mousePointer);
    }
}

void SoundDrawWindow::mouseMoved(QMouseEvent* e)
{
    if (nullptr == mSelectionRect) return;
    const int x = e->x(), y = e->y();
    const QRect CR = Bild->contentsRect();
    const double left = mSelectionRect->left() * XSF + XO, top = YO - mSelectionRect->top() * YSF, right = mSelectionRect->right() * XSF + XO, bottom = YO - mSelectionRect->bottom() * YSF;
    if (mMoveState == MSOutside)
    {
        if (x < left - 5 || x < CR.left() + ScaleYWidth || y < top - 5 || y < CR.top() || x > right + 5 || x > CR.right() || y > bottom + 5 || y > CR.bottom() - ScaleXHeight)
        {
            ensureMouseShape(Qt::ArrowCursor);
            mMouseState = MSOutside;
        }
        else
        {
            if (x > left + 5 && y > top + 5 && x < right - 5 && y < bottom - 5)
            {
                ensureMouseShape(Qt::SizeAllCursor);
                mMouseState = MSInside;
            }
            else if (x <= left + 5)
            {
                if (y <= top + 5)
                {
                    ensureMouseShape(Qt::SizeFDiagCursor);
                    mMouseState = MSLTCorner;
                }
                else if (y >= bottom - 5)
                {
                    ensureMouseShape(Qt::SizeBDiagCursor);
                    mMouseState = MSBLCorner;
                }
                else
                {
                    ensureMouseShape(Qt::SizeHorCursor);
                    mMouseState = MSLeft;
                }
            }
            else if (x >= right - 5)
            {
                if (y <= top + 5)
                {
                    ensureMouseShape(Qt::SizeBDiagCursor);
                    mMouseState = MSTRCorner;
                }
                else if (y >= bottom - 5)
                {
                    ensureMouseShape(Qt::SizeFDiagCursor);
                    mMouseState = MSRBCorner;
                }
                else
                {
                    ensureMouseShape(Qt::SizeHorCursor);
                    mMouseState = MSRight;
                }
            }
            else
            {
                ensureMouseShape(Qt::SizeVerCursor);
                mMouseState = (y >= bottom - 5 ? MSBottom : MSTop);
            }
            e->accept();
        }
    }
    else
    {
        const double mouseXDiff = x - mMoveMouseStartPoint.x(), mouseYDiff = y - mMoveMouseStartPoint.y();
        switch (mMouseState)
        {
            case MSInside:
                mSelectionRect->setLeft(mMoveStartRect.left() + mouseXDiff / XSF);
                mSelectionRect->setTop(mMoveStartRect.top() - mouseYDiff / YSF);
                mSelectionRect->setWidth(mMoveStartRect.width());
                mSelectionRect->setHeight(mMoveStartRect.height());
                break;
            case MSLeft:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                    mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                }
                break;
            case MSTop:
                {
                    const double yMoveObcC = mouseYDiff / YSF;
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSLTCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    const double yMoveObcC = mouseYDiff / YSF;
                    mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                    mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSTRCorner:
                {
                    const double yMoveObcC = mouseYDiff / YSF;
                    mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSRight:
                mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                break;
            case MSRBCorner:
                mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                mSelectionRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBottom:
                mSelectionRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBLCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                    mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    mSelectionRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                }
                break;
            case MSOutside:
            default:
                // don't get here
                break;
        }
        e->accept();
        Paint();
    }
}

void SoundDrawWindow::mousePressed(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton)
    {
        mMoveState = mMouseState;
        if (mMouseState != MSOutside)
        {
            mMoveMouseStartPoint = QPoint(e->x(), e->y());
            mMoveStartRect = *mSelectionRect;
        }
        ensureMouseShape(Qt::DragMoveCursor);
    }
    else if (e->button() == Qt::RightButton) ShowPopupMenu(e->globalPos());
}

void SoundDrawWindow::mouseReleased(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton)
    {
        mMoveState = MSOutside;
        ensureMouseShape(Qt::ArrowCursor);
    }
}

void SoundDrawWindow::ShowPopupMenu(const QPoint& point)
{
    mPopupMenu->popup(point);
}

int SoundDrawWindow::getSoundData(float ** data)
{
    int xStart = 0, length = Daten->GetDSL(), xStop = length - 1, n = 1;
    if (nullptr != mSelectionRect)
    {
        const double startTime = mSelectionRect->left(), stopTime = mSelectionRect->right();
        while (n < length && Daten->GetValue(n, 0) < startTime) ++n;
        xStart = n;
        while (n < length && Daten->GetValue(n, 0) <= stopTime) ++n;
        xStop = n-1;
    }
    int rLength = xStop - xStart + 1, i;
    *data = new float[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = static_cast<float>(Daten->GetValue(n, 1));
    return rLength;
}

void SoundDrawWindow::Play()
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
}

void SoundDrawWindow::WriteToFile()
{
    QString filename = QFileDialog::getSaveFileName(this, "Select filename to save");
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    float* data;
    int length = getSoundData(&data);
    file.write(reinterpret_cast<char*>(data), 4*length);
}
