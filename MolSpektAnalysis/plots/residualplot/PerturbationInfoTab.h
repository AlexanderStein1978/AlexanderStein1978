//
// C++ Interface: PerturbationInfoTab
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef PERTURBATIONINFOTAB_H
#define PERTURBATIONINFOTAB_H


#include <QTableWidget>


class PerturbationInfoTab : public QTableWidget
{
    Q_OBJECT

public:

    PerturbationInfoTab(QWidget* parent = 0);

protected:

    virtual void mouseReleaseEvent(QMouseEvent *event);

signals:

    void RightClicked(QPoint pos);
};

#endif
