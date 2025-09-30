//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "PerturbationInfoTab.h"

#include <QMouseEvent>


PerturbationInfoTab::PerturbationInfoTab(QWidget *parent) : QTableWidget(parent)
{
}

void PerturbationInfoTab::mouseReleaseEvent(QMouseEvent *event)
{
    QTableWidget::mouseReleaseEvent(event);
    if (event->button() == Qt::RightButton)
    {
        event->accept();
        emit RightClicked(mapToGlobal(event->pos()));
    }
}
