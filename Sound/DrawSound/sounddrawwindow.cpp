#include "sounddrawwindow.h"

#include "datensatz.h"

#include <QPainter>
#include <QToolButton>
#include <QComboBox>
#include <QMenu>
#include <QFileDialog>

#include <cmath>


SoundDrawWindow::SoundDrawWindow(SoundRecordAndDrawControl *const control, const int sampleRate, const int o) : DiagWindow(SimpleDiagWindow, nullptr, "Data files (*.dat)", ".dat", o),
    mPopupMenu(new QMenu(this)), mSampleRate(sampleRate), mControl(control)
{
    QAction *writeAct = new QAction("Write to file...", this);
    mPopupMenu->addAction(writeAct);
    connect(writeAct, SIGNAL(triggered()), this, SLOT(WriteToFile()));
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(SelectionChanged(QRect*)));
    connect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
    connect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
    connect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
}

SoundDrawWindow::~SoundDrawWindow()
{
    if (nullptr != mSelectionRect) delete mSelectionRect;
}

void SoundDrawWindow::SelectionChanged(QRect* MarkedArea)
{
    if (!ZoomB->isChecked())
    {
        if (nullptr == mSelectionRect) mSelectionRect = new QRectF;
        if (mIsFFT)
        {
            double left = double(MarkedArea->left() - XO) / XSF, right = double(MarkedArea->right() - XO) / XSF, width = right - left, FFTWidth = getFFTWidth(width), delta = width - FFTWidth;
            mSelectionRect->setCoords(left + delta, -double(MarkedArea->top() - YO) / YSF, right - delta, -double(MarkedArea->bottom() - YO) / YSF);
            showFFT();
        }
        else mSelectionRect->setCoords(double(MarkedArea->left() - XO) / XSF, -double(MarkedArea->top() - YO) / YSF, double(MarkedArea->right() - XO) / XSF, -double(MarkedArea->bottom() - YO) / YSF);
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
    if (nullptr == mSelectionRect || ZoomB->isChecked()) return;
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
                if (mIsFFT) showFFT();
                break;
            case MSLeft:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mSelectionRect->setLeft(mSelectionRect->right() - newWidth);
                        mSelectionRect->setWidth(newWidth);
                        if (oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
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
                    if (mIsFFT)
                    {
                        const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mSelectionRect->setLeft(mSelectionRect->right() - newWidth);
                        mSelectionRect->setWidth(newWidth);
                        if (oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSTRCorner:
                {
                    const double yMoveObcC = mouseYDiff / YSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                        mSelectionRect->setWidth(newWidth);
                        if (oldWidth != newWidth) showFFT();
                    }
                    else mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSRight:
                if (mIsFFT)
                {
                    const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mSelectionRect->setWidth(newWidth);
                    if (oldWidth != newWidth) showFFT();
                }
                else mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                break;
            case MSRBCorner:
                if (mIsFFT)
                {
                    const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mSelectionRect->setWidth(newWidth);
                    if (oldWidth != newWidth) showFFT();
                }
                else mSelectionRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                mSelectionRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBottom:
                mSelectionRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBLCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mSelectionRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mSelectionRect->setLeft(mSelectionRect->right() - newWidth);
                        mSelectionRect->setWidth(newWidth);
                        if (oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
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
    if (!ZoomB->isChecked() &&  e->button() == Qt::LeftButton)
    {
        mMoveState = mMouseState;
        if (mMouseState != MSOutside)
        {
            mMoveMouseStartPoint = QPoint(e->x(), e->y());
            mMoveStartRect = *mSelectionRect;
        }
        ensureMouseShape(Qt::DragMoveCursor);
    }
}

void SoundDrawWindow::mouseReleased(QMouseEvent* e)
{
    if (!ZoomB->isChecked())
    {
        if (e->button() == Qt::LeftButton)
        {
            mMoveState = MSOutside;
            ensureMouseShape(Qt::ArrowCursor);
        }
        else if (e->button() == Qt::RightButton) ShowPopupMenu(e->globalPos());
    }
}

void SoundDrawWindow::ShowPopupMenu(const QPoint& point)
{
    mPopupMenu->popup(point);
}

int SoundDrawWindow::getSoundDataRange(int& xStart, int& xStop)
{
    int length = Daten->GetDSL(), n = 1;
    xStart = 0;
    xStop = length - 1;
    if (nullptr != mSelectionRect)
    {
        const double startTime = mSelectionRect->left(), stopTime = mSelectionRect->right();
        while (n < length && Daten->GetValue(n, 0) < startTime) ++n;
        xStart = n;
        while (n < length && Daten->GetValue(n, 0) <= stopTime) ++n;
        xStop = n-1;
    }
    return xStop - xStart + 1;
}

double SoundDrawWindow::getFFTWidth(const double inputWidth)
{
    int N = static_cast<int>(log2(inputWidth * mSampleRate));
    if (N < 1) N=1;
    double rc = (static_cast<double>(2 << N)) / mSampleRate;
    return (abs(inputWidth - rc) < abs(inputWidth - 0.5 * rc) ? rc : 0.5 * rc);
}
