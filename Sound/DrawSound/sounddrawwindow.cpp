#include "sounddrawwindow.h"

#include <QPainter>
#include <QToolButton>


SoundDrawWindow::SoundDrawWindow(QString filename)
{
    setUnits("time [s]", "intensity");
    setWindowTitle("Draw sound: " + filename);
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(SelectionChanged(QRect*)));
    connect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
    connect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
    connect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
}

SoundDrawWindow::~SoundDrawWindow() noexcept
{
    if (nullptr != mSelectionRect) delete mSelectionRect;
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
                    mMouseState = MSBLCorner;
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
                    mSelectionRect->setWidth(mMoveStartRect.width() - mouseXDiff / XSF);
                    mSelectionRect->setTop(mMoveStartRect.top() - yMoveObcC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + yMoveObcC);
                }
                break;
            case MSRight:
                mSelectionRect->setWidth(mMoveStartRect.width() - mouseXDiff / XSF);
                break;
            case MSRBCorner:
                mSelectionRect->setWidth(mMoveStartRect.width() - mouseXDiff / XSF);
                mSelectionRect->setHeight(mMoveStartRect.height() + mouseYDiff / YSF);
                break;
            case MSBottom:
                mSelectionRect->setHeight(mMoveStartRect.height() + mouseYDiff / YSF);
                break;
            case MSBLCorner:
                {
                    const double xMoveObjC = mouseXDiff / XSF;
                    mSelectionRect->setLeft(mMoveStartRect.left() + xMoveObjC);
                    mSelectionRect->setWidth(mMoveStartRect.width() - xMoveObjC);
                    mSelectionRect->setHeight(mMoveStartRect.height() + mouseYDiff / YSF);
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
}

void SoundDrawWindow::mouseReleased(QMouseEvent* e)
{
    mMoveState = MSOutside;
    ensureMouseShape(Qt::ArrowCursor);
}
