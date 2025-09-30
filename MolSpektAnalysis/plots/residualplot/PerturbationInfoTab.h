//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
