//
// C++ Implementation: TSDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "tsdialog.h"
#include "molecule.h"
#include "linetable.h"

#include <QListWidget>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>


TSDialog::TSDialog(QWidget* parent, Molecule* M, ElState* S)
    : QDialog(parent)
    , Tables(0)
    , NTables(0)
    , Mol(M)
{
    int i, n, N = M->getNumTransitions();
    Transition *T = 0;
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Electronic transition:", this), 0, 0);
    L->addWidget(TransBox = new QComboBox(this), 0, 1);
    QLabel *L1 = new QLabel("Selected tables:", this);
    L->addWidget(L1, 1, 0);
    QLabel *L2 = new QLabel("Available tables:", this);
    L->addWidget(L2, 1, 1);
    List = new QListWidget(this);
    L->addWidget(List, 2, 0, 3, 1);
    Box = new QComboBox(this);
    L->addWidget(Box, 2, 1);
    Add = new QPushButton("Add table", this);
    L->addWidget(Add, 3, 1);
    Remove = new QPushButton("Remove table", this);
    L->addWidget(Remove, 4, 1);
    QGridLayout *La2 = new QGridLayout;
    KeepOther = new QCheckBox("Keep other data", this);
    KeepOther->setChecked(true);
    La2->addWidget(KeepOther, 0, 0);
    UpdateSA = new QCheckBox("Update state assignment", this);
    La2->addWidget(UpdateSA, 0, 1);
    L->addLayout(La2, 5, 0, 2, 1);
    L->setRowMinimumHeight(6, 20);
    OKB = new QPushButton("OK", this);
    L->addWidget(OKB, 7, 0);
    Cancel = new QPushButton("Cancel", this);
    L->addWidget(Cancel, 7, 1);
    for (n=0; n<N; ++n) TransBox->addItem(M->getTransitionP(n)->getName());
    for (n=0; n<N; ++n) if ((T=M->getTransitionP(n))->getLowerState() == S || T->getUpperState() == S)
    {
        TransBox->setCurrentIndex(n);
        break;
    }
    TransBox->setEditable(false);
    if (n<N) for (i=0; i < T->getNumLineTables(); i++) Box->addItem(T->getLineTableName());
    Box->setEditable(false);
    connect(TransBox, SIGNAL(currentIndexChanged(int)), this, SLOT(transBoxChanged(int)));
    connect(KeepOther, SIGNAL(toggled(bool)), this, SLOT(KeepOtherChanged(bool)));
    connect(Add, SIGNAL(clicked()), this, SLOT(addTable()));
    connect(Remove, SIGNAL(clicked()), this, SLOT(removeTable()));
    connect(OKB, SIGNAL(clicked()), this, SLOT(OK()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

TSDialog::~TSDialog()
{
}

void TSDialog::addTable()
{
    int n;
    if ((n = Box->currentIndex()) != -1 && List->findItems(Box->currentText(), Qt::MatchExactly).isEmpty()) List->addItem(Box->currentText());
}

bool TSDialog::getKeepOther()
{
    return KeepOther->isChecked();
}

int TSDialog::getNumSelected()
{
    return NTables;
}

LineTable* TSDialog::getSelected(int n)
{
    return Tables[n];
}

void TSDialog::OK()
{
    int n, m, N = List->count();
    bool S = true;
    Tables = new LineTable*[N];
    for (n=m=0; n<N; n++)
    {
        Tables[m] = Mol->getLineTable(List->item(m)->text());
        if (Tables[m] != 0) ++m;
        else
        {
            S = false;
            delete List->takeItem(m);
        }
    }
    NTables = m;
    if(S) accept();
}

void TSDialog::removeTable()
{
    int n;
    if ((n = List->currentRow()) != -1) delete List->takeItem(n);
}

void TSDialog::transBoxChanged(int i)
{
    Box->clear();
    Transition *T = Mol->getTransitionP(i);
    int n, N = T->getNumLineTables();
    for (n=0; n<N; ++n) Box->addItem(T->getLineTableName(n));
}
