#pragma once

#include "DiagWindow.h"


class SoundDrawWindow : public DiagWindow
{
    Q_OBJECT
public:
    SoundDrawWindow(QString filename);
    ~SoundDrawWindow();

private slots:
    void SelectionChanged(QRect *MarkedArea);
    void mouseMoved(QMouseEvent *e);
    void mousePressed(QMouseEvent *e);
	void mouseReleased(QMouseEvent *e);

private:
    enum MouseState{MSOutside, MSInside, MSLeft, MSLTCorner, MSTop, MSTRCorner, MSRight, MSRBCorner, MSBottom, MSBLCorner};

    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    void ensureMouseShape(const Qt::CursorShape shape);

    QRectF *mSelectionRect = nullptr, mMoveStartRect;
    MouseState mMouseState = MSOutside, mMoveState = MSOutside;
    QPoint mMoveMouseStartPoint;
};
