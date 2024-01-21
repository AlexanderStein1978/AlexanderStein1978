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

private:
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    void ensureMouseShape(const Qt::CursorShape shape);

    QRectF* mSelectionRect = nullptr;
};
