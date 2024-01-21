#include "sounddrawwindow.h"

#include <QPainter>
#include <QToolButton>


SoundDrawWindow::SoundDrawWindow(QString filename)
{
    setUnits("time [s]", "intensity");
    setWindowTitle("Draw sound: " + filename);
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(SelectionChanged(QRect*)));
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
    if (x < left - 5 || x < CR.left() + ScaleYWidth || y < top - 5 || y < CR.top() || x > right + 5 || x > CR.right() || y > bottom + 5 || x > CR.bottom() - ScaleXHeight)
        ensureMouseShape(Qt::ArrowCursor);
    else
    {
        if (x > left + 5 && y > top + 5 && x < right - 5 && y < bottom - 5) ensureMouseShape(Qt::SizeAllCursor);
        else if (x <= left + 5)
        {
            if (y <= top + 5) ensureMouseShape(Qt::SizeFDiagCursor);
            else if (y >= bottom - 5) ensureMouseShape(Qt::SizeBDiagCursor);
            else ensureMouseShape(Qt::SizeHorCursor);
        }
        else if (x >= right - 5)
        {
            if (y <= top + 5) ensureMouseShape(Qt::SizeBDiagCursor);
            else if (y >= bottom - 5) ensureMouseShape(Qt::SizeFDiagCursor);
            else ensureMouseShape(Qt::SizeHorCursor);
        }
        else ensureMouseShape(Qt::SizeVerCursor);
        e->accept();
    }
}
