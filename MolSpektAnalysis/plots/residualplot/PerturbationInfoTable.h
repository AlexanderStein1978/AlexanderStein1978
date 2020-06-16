//
// C++ Interface: PerturbationInfoTable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef PERTURBATIONINFOTABLE_H
#define PERTURBATIONINFOTABLE_H


#include <QWidget>

class ResidualFit;
class PerturbationInfoTab;

class QComboBox;


class PerturbationInfoTable : public QWidget
{
    Q_OBJECT

public:

    PerturbationInfoTable(QWidget* i_parent = 0);
    ~PerturbationInfoTable();

    void SetResidualFit(ResidualFit* const i_ResidualFit);

protected:

    virtual void mouseReleaseEvent(QMouseEvent *event);

public slots:

    void SelectPerturbation(int i_index);

private slots:

    void ReduceNumberOfBandConstants();
    void ShowPopupMenu(QPoint i_pos);

private:

    PerturbationInfoTab* m_Tab;
    QComboBox* m_PerturbationSelectBox;
    ResidualFit* m_ResidualFit;
};

#endif
