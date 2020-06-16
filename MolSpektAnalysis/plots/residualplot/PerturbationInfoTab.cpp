//
// C++ Implementation: PerturbationInfoTab
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
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
