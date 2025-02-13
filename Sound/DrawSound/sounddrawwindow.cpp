#include "sounddrawwindow.h"
#include "soundmainwindow.h"

#include "datensatz.h"

#include <QPainter>
#include <QToolButton>
#include <QComboBox>
#include <QMenu>
#include <QFileDialog>
#include <QLineEdit>

#include <cmath>


SoundDrawWindow::SoundDrawWindow(SoundRecordAndDrawControl *const control, SoundMainWindow *const MW, const int sampleRate, const int o) : DiagWindow(SimpleDiagWindow, MW, "Data files (*.dat)", ".dat", o),
    mPopupMenu(new QMenu(this)), mSampleRate(sampleRate), mControl(control), mMode(MNormal)
{
    QAction *writeAct = new QAction("Write to file...", this);
    mPopupMenu->addAction(writeAct);
    connect(writeAct, SIGNAL(triggered()), this, SLOT(WriteToFile()));
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(SelectionChanged(QRect*)));
    connect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
    connect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
    connect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
    connect(Bild, SIGNAL(LeftClicked(QPoint*)), this, SLOT(mouseLeftClicked(QPoint*)));
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
            double left = double(MarkedArea->left() - XO) / XSF, right = double(MarkedArea->right() - XO) / XSF, width = right - left, FFTWidth = getFFTWidth(width), delta = 0.5 * (width - FFTWidth);
            mSelectionRect->setCoords(left + delta, -double(MarkedArea->top() - YO) / YSF, right - delta, -double(MarkedArea->bottom() - YO) / YSF);
            showFFT();
        }
        else mSelectionRect->setCoords(double(MarkedArea->left() - XO) / XSF, -double(MarkedArea->top() - YO) / YSF, double(MarkedArea->right() - XO) / XSF, -double(MarkedArea->bottom() - YO) / YSF);
        updateLabelRectSelections();
        if (mMode == MFastLabeling && ((mSelectionRect->left() <= mXStart && mXStart > 0 && mSelectionRect->right() < mXStop) ||
                                       (mSelectionRect->right() >= mXStop && mXStop < Daten->GetValue(Daten->GetDSL() - 1, 0) && mSelectionRect->left() > mXStop)))
        {
            double diff = 0.5 * (mSelectionRect->left() + mSelectionRect->right() - mXStart - mXStop);
            if (mSelectionRect->left() + diff < 0) diff = -1.0 * mSelectionRect->left();
            else if (mSelectionRect->right() + diff > Daten->GetValue(Daten->GetDSL() - 1, 0)) diff = Daten->GetValue(Daten->GetDSL() - 1, 0) - mSelectionRect->right();
            xStart->setText(QString::number(mXStart + diff, 'g', 11));
            xStop->setText(QString::number(mXStop + diff, 'g', 11));
            QPoint mousePos = this->mapFromGlobal(QCursor::pos());
            mousePos.setX(mousePos.x()  - diff * XSF);
            QCursor::setPos(this->mapToGlobal(mousePos));
        }
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

void SoundDrawWindow::mouseLeftClicked(QPoint* point)
{
    std::pair<int, MouseState> state = getBestMouseState(*point);
    if ((state.second == MSOutside || state.first != -1) && nullptr != mSelectionRect)
    {
        delete mSelectionRect;
        mSelectionRect = nullptr;
    }
    if (state.second != MSOutside) for (int i=0; i < mLabels.size(); ++i) mLabels[i].isSelected = (i == state.first);
    else for (int i=0; i < mLabels.size(); ++i) mLabels[i].isSelected = false;
    Paint();
}

void SoundDrawWindow::mouseMoved(QMouseEvent* e)
{
    if ((nullptr == mSelectionRect && 0 == mLabels.size()) || ZoomB->isChecked()) return;
    if (mMoveState == MSOutside)
    {
        mMouseState = getBestMouseState(e->pos()).second;
        switch (mMouseState)
        {
            case MSBLCorner:
            case MSTRCorner:
               ensureMouseShape(Qt::SizeBDiagCursor);
               break;
            case MSBottom:
            case MSTop:
                ensureMouseShape(Qt::SizeVerCursor);
                break;
            case MSInside:
                ensureMouseShape(Qt::SizeAllCursor);
                break;
            case MSLeft:
            case MSRight:
                ensureMouseShape(Qt::SizeHorCursor);
                break;
            case MSLTCorner:
            case MSRBCorner:
                ensureMouseShape(Qt::SizeFDiagCursor);
                break;
            case MSOutside:
                ensureMouseShape(Qt::ArrowCursor);
                break;
        }
        if (mMouseState != MSOutside) e->accept();
    }
    else
    {
        const int x = e->x(), y = e->y();
        const double mouseXDiff = x - mMoveMouseStartPoint.x(), mouseYDiff = y - mMoveMouseStartPoint.y();
        switch (mMouseState)
        {
            case MSInside:
                mMovingRect->setLeft(mMoveStartRect.left() + mouseXDiff / XSF);
                mMovingRect->setTop(mMoveStartRect.top() - mouseYDiff / YSF);
                mMovingRect->setWidth(mMoveStartRect.width());
                mMovingRect->setHeight(mMoveStartRect.height());
                if (mIsFFT && mMovingRect == mSelectionRect) showFFT();
                break;
            case MSLeft:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mMovingRect->setLeft(mMovingRect->right() - newWidth);
                        mMovingRect->setWidth(newWidth);
                        if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mMovingRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mMovingRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
                }
                break;
            case MSTop:
                {
                    const double yMoveObcC = mouseYDiff / YSF;
                    mMovingRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mMovingRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSLTCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    const double yMoveObcC = mouseYDiff / YSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mMovingRect->setLeft(mMovingRect->right() - newWidth);
                        mMovingRect->setWidth(newWidth);
                        if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mMovingRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mMovingRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
                    mMovingRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mMovingRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSTRCorner:
                {
                    const double yMoveObcC = mouseYDiff / YSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                        mMovingRect->setWidth(newWidth);
                        if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                    }
                    else mMovingRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mMovingRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mMovingRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSRight:
                if (mIsFFT)
                {
                    const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mMovingRect->setWidth(newWidth);
                    if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                }
                else mMovingRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                break;
            case MSRBCorner:
                if (mIsFFT)
                {
                    const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                    mMovingRect->setWidth(newWidth);
                    if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                }
                else mMovingRect->setWidth(mMoveStartRect.width() + mouseXDiff / XSF);
                mMovingRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBottom:
                mMovingRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                break;
            case MSBLCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    if (mIsFFT)
                    {
                        const double oldWidth = mMovingRect->width(), newWidth = getFFTWidth(mMoveStartRect.width() - xMoveObjC);
                        mMovingRect->setLeft(mMovingRect->right() - newWidth);
                        mMovingRect->setWidth(newWidth);
                        if (mMovingRect == mSelectionRect && oldWidth != newWidth) showFFT();
                    }
                    else
                    {
                        mMovingRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                        mMovingRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    }
                    mMovingRect->setHeight(mMoveStartRect.height() - mouseYDiff / YSF);
                }
                break;
            case MSOutside:
            default:
                // don't get here
                break;
        }
        e->accept();
        if (mMovingRect == mSelectionRect) updateLabelRectSelections();
        Paint();
    }
}

void SoundDrawWindow::mousePressed(QMouseEvent* e)
{
    if (!ZoomB->isChecked() &&  e->button() == Qt::LeftButton)
    {
        std::pair<int, MouseState> state = getBestMouseState(e->pos());
        mMoveState = mMouseState = state.second;
        if (mMouseState != MSOutside)
        {
            mMovingRect = (state.first >= 0 ? &mLabels[state.first].rect : mSelectionRect);
            mMoveMouseStartPoint = QPoint(e->x(), e->y());
            mMoveStartRect = *mMovingRect;
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

int SoundDrawWindow::getSoundDataRange(int& xStart, int& xStop, const int labelIndex, const FFTSelection fftSelection) const
{
    const QRectF* rect = (labelIndex >= 0 ? &mLabels[labelIndex].rect : mSelectionRect);
    double start =  (nullptr != rect ? rect->left() : 0), stop = (nullptr != rect ? rect->right() : Daten->GetValue(Daten->GetDSL() - 1, 0));
    if (labelIndex >= 0 && (fftSelection == FSForFTT  || fftSelection == FSForSelectedFTT || (fftSelection == FSDependsOnState && mIsFFT)))
    {
        double length = stop - start, fftLength = (fftSelection == FSForSelectedFTT ? mSelectedFFtSize : getFFTWidth(length));
        if (fftLength > length) fftLength *= 0.5;
        double step = 0.5 * (length - fftLength);
        start += step;
        stop -= step;
    }
    int length = Daten->GetDSL(), n = 1;
    while (n < length && Daten->GetValue(n, 0) < start) ++n;
    xStart = n;
    while (n < length && Daten->GetValue(n, 0) <= stop) ++n;
    xStop = n-1;
    return xStop - xStart + 1;
}

int SoundDrawWindow::getFFTLength(const int inputLength) const
{
    int N = static_cast<int>(log2(inputLength));
    if (N < 1) N=1;
    int rc = (2 << N);
    return (abs(inputLength - rc) < abs(inputLength - rc / 2) ? rc : rc / 2);
}

double SoundDrawWindow::getFFTWidth(const double inputWidth) const
{
    return static_cast<double>(getFFTLength(static_cast<int>(inputWidth * mSampleRate))) / mSampleRate;
}

void SoundDrawWindow::updateLabelRectSelections()
{
    const double left = mSelectionRect->left() * XSF + XO, top = YO - mSelectionRect->top() * YSF, right = mSelectionRect->right() * XSF + XO, bottom = YO - mSelectionRect->bottom() * YSF;
    for (Label label : mLabels) label.isSelected = (left <= label.rect.left() && right >= label.rect.right() && top <= label.rect.top() && bottom >= label.rect.top());
}

std::pair<int, SoundDrawWindow::MouseState> SoundDrawWindow::getBestMouseState(const QPoint& point)
{
    std::pair<int, MouseState> rValue;
    rValue.first = -1;
    for (int i=0; i < mLabels.size(); ++i) if (getLabelTextRect(mLabels[i]).contains(point))
    {
        rValue.first = i;
        rValue.second = MSInside;
        return rValue;
    }
    if (nullptr != mSelectionRect && (rValue.second = isCloseToCorner(*mSelectionRect, point)) != MSOutside) return rValue;
    for (int i=0; i < mLabels.size(); ++i) if ((rValue.second = isCloseToCorner(mLabels[i].rect, point)) != MSOutside)
    {
        rValue.first = i;
        return rValue;
    }
    if (nullptr != mSelectionRect && (rValue.second = isCloseToWall(*mSelectionRect, point)) != MSOutside) return rValue;
    for (int i=0; i < mLabels.size(); ++i) if ((rValue.second = isCloseToWall(mLabels[i].rect, point)) != MSOutside)
    {
        rValue.first = i;
        return rValue;
    }
    if (nullptr != mSelectionRect && isInsideRect(*mSelectionRect, point)) rValue.second = MSInside;
    return rValue;
}

QRect SoundDrawWindow::getLabelTextRect(const Label& label)
{
    int bottom = YO - label.rect.bottom() * YSF, height = label.rect.height() * YSF, left = label.rect.left() * XSF + XO, width = label.rect.width() * XSF;
    int rWidth = TextWidth(mLabelFont, label.phoneme), rHeight = TextHeight(mLabelFont, label.phoneme);
    return QRect(left + (width - rWidth) / 2, bottom + height - (5 * rHeight) / 4, rWidth, rHeight);
}

bool SoundDrawWindow::arePointsClose(const QPointF& pointF, const QPoint& point) const
{
    const int left = static_cast<int>(pointF.x() * XSF + XO), top = static_cast<int>(YO - pointF.y() * YSF);
    return abs(left - point.x()) <= D && abs(top - point.y()) <= D;
}

int SoundDrawWindow::estimateLabelIndex(const QString& phoneme)
{
    int index = mLabelOrder.indexOf(phoneme);
    if (-1 == index)
    {
        index = mLabelOrder.size();
        mLabelOrder.push_back(phoneme);
    }
    return index;
}

SoundDrawWindow::MouseState SoundDrawWindow::isCloseToCorner(const QRectF& rect, const QPoint& point) const
{
    if (arePointsClose(rect.bottomLeft(), point)) return MSBLCorner;
    if (arePointsClose(rect.bottomRight(), point)) return MSRBCorner;
    if (arePointsClose(rect.topLeft(), point)) return MSLTCorner;
    if (arePointsClose(rect.topRight(), point)) return MSTRCorner;
    return MSOutside;
}

SoundDrawWindow::MouseState SoundDrawWindow::isCloseToWall(const QRectF& rect, const QPoint& point) const
{
    const int left = static_cast<int>(rect.left() * XSF + XO), top = static_cast<int>(YO - rect.top() * YSF), right = static_cast<int>(rect.right() * XSF + XO);
    const int bottom = static_cast<int>(YO - rect.bottom() * YSF), x = point.x(), y = point.y();
    if (y > top && y < bottom)
    {
        if (abs(x - left) <= D) return MSLeft;
        if (abs(x - right) <= D) return MSRight;
    }
    if (x > left && x < right)
    {
        if (abs(y - top) <= D) return MSTop;
        if (abs(y - bottom) <= D) return MSBottom;
    }
    return MSOutside;
}

bool SoundDrawWindow::isInsideRect(const QRectF& rect, const QPoint& point) const
{
    const int left = static_cast<int>(rect.left() * XSF + XO), top = static_cast<int>(YO - rect.top() * YSF), right = static_cast<int>(rect.right() * XSF + XO);
    const int bottom = static_cast<int>(YO - rect.bottom() * YSF), x = point.x(), y = point.y();
    return (x >= left && x <= right && y >= top && y <= bottom);
}
