//
// C++ Implementation: EsPotFitInputElStateAssignDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "espotfitinputelstateassigndialog.h"
#include "elstate.h"

#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QPushButton>


EsPotFitInputElStateAssignDialog::EsPotFitInputElStateAssignDialog(QWidget* parent, ElState** i_states, int i_numStates)
    : QDialog(parent), m_boxes(new QComboBox*[i_numStates])
    , m_assignmentFromTab(new QCheckBox("Use assignment from vss column", this)), m_numStates(i_numStates)
{

    QGridLayout *L = new QGridLayout(this);
    setWindowTitle("Please select electronic state assignment");
    L->addWidget(m_assignmentFromTab, 0, 0, 1, 2);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Electronic State:", this), 2, 0);
    L->addWidget(new QLabel("Assignment:", this), 2, 1);
    L->setRowMinimumHeight(3, 10);
    for (int i=0; i < i_numStates; ++i)
    {
        L->addWidget(new QLabel(i_states[i]->getName(), this), i+4, 0);
        L->addWidget(m_boxes[i] = new QComboBox(this), i+4, 1);
        for (int s=i; s < i_numStates; ++s) m_boxes[i]->addItem(QString::number(-1 - s));
        m_boxes[i]->setEditable(false);
        connect(m_boxes[i], SIGNAL(currentIndexChanged(int)), this, SLOT(BoxSelectionChanged()));
    }
    L->setRowMinimumHeight(i_numStates + 4, 20);
    QPushButton *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    L->addWidget(OK, i_numStates + 5, 0);
    L->addWidget(Cancel, i_numStates + 5, 1);
    connect(m_assignmentFromTab, SIGNAL(toggled(bool)), this, SLOT(AssignmentFromTabSelected(bool)));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

EsPotFitInputElStateAssignDialog::~EsPotFitInputElStateAssignDialog()
{
    delete[] m_boxes;
}

int* EsPotFitInputElStateAssignDialog::GetAssignment() const
{
    int *returnArray = new int[m_numStates];
    for (int i=0; i < m_numStates; ++i) returnArray[i] = m_boxes[i]->currentText().toInt();
    return returnArray;
}

void EsPotFitInputElStateAssignDialog::AssignmentFromTabSelected(bool i_selected)
{
    for (int i=0; i < m_numStates; ++i) m_boxes[i]->setEnabled(!i_selected);
}

void EsPotFitInputElStateAssignDialog::BoxSelectionChanged()
{
    bool AssignmentsToUse[m_numStates];
    int s, i;
    for (i=0; i < m_numStates; ++i) AssignmentsToUse[i] = true;
    for (i=1; i < m_numStates; ++i)
    {
        AssignmentsToUse[-1 - m_boxes[i-1]->currentText().toInt()] = false;
        if (!AssignmentsToUse[-1 - m_boxes[i]->currentText().toInt()])
        {
            m_boxes[i]->clear();
            for (s=0; s < m_numStates; ++s) if (AssignmentsToUse[s]) m_boxes[i]->addItem(QString::number(-1 - s));
        }
    }
}
