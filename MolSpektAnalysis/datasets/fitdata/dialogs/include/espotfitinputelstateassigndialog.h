//
// C++ Interface: EsPotFitInputElStateAssignDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef ESPOTFITINPUTELSTATEASSIGNDIALOG_H
#define ESPOTFITINPUTELSTATEASSIGNDIALOG_H


#include <QDialog>
#include <QCheckBox>

class QComboBox;

class ElState;


class EsPotFitInputElStateAssignDialog : public QDialog
{
    Q_OBJECT

public:
    EsPotFitInputElStateAssignDialog(QWidget* parent, ElState** i_StateArray, int i_numStates);
    virtual ~EsPotFitInputElStateAssignDialog();
    int* GetAssignment() const;

    inline bool IsAssignmentFromTabSelected() const
    {
        return m_assignmentFromTab->isChecked();
    }

private slots:
    void AssignmentFromTabSelected(bool i_selected);
    void BoxSelectionChanged();

private:
    QComboBox** m_boxes;
    QCheckBox* m_assignmentFromTab;
    int m_numStates;
};

#endif
