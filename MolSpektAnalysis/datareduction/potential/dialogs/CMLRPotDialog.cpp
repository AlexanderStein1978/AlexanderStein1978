//
// C++ Implementation: CMLRPotDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "CMLRPotDialog.h"


CMLRPotDialog::CMLRPotDialog(QStringList *nCVL, QWidget* parent): QDialog(parent)
{
    CVL = nCVL;
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("p", this), 0, 0);
    L->addWidget(pE = new QLineEdit("5", this), 0, 1);
    L->addWidget(new QLabel("q", this), 0, 2);
    L->addWidget(qE = new QLineEdit("3", this), 0, 3);
    L->addWidget(new QLabel("R_ref [A]", this), 0, 4);
    L->addWidget(RrefE = new QLineEdit("5.5", this), 0, 5);
    L->addWidget(new QLabel("N beta coeff.:", this), 0, 6);
    L->addWidget(NCE = new QLineEdit("10", this), 0, 7);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Long range coefficients:", this), 2, 0, 1, 8);
    L->addWidget(LRC = new QListWidget(this), 3, 0, 3, 4);
    L->addWidget(addB = new QPushButton("<- add coefficient", this), 3, 4, 1, 2);
    L->addWidget(new QLabel("new LRC:", this), 3, 6);
    L->addWidget(LRCB = new QComboBox(this), 3, 7);
    LRCB->addItems(QStringList() << "C3" << "C6" << "C8" << "C10" << "C12");
    LRCB->setEditable(false);
    L->addWidget(delB = new QPushButton("<- delete coeff.", this), 4, 4, 1, 2);
    L->addWidget(setB = new QPushButton("<- set value", this), 5, 4, 1, 2);
    L->addWidget(new QLabel("Value [cm^-1/A^n]:", this), 5, 6);
    L->addWidget(CValE = new QLineEdit("1.525e7", this), 5, 7);
    L->setRowMinimumHeight(6, 20);
    L->addWidget(new QLabel("Fit range:", this), 7, 0, 1, 8);
    QGridLayout *L2 = new QGridLayout();
    L2->addWidget(new QLabel("R_min [A]:", this), 0, 0);
    L2->addWidget(mRE = new QLineEdit("3.0", this), 0, 1);
    L2->addWidget(new QLabel("R_max [A]:", this), 0, 2);
    L2->addWidget(MRE = new QLineEdit("30,0", this), 0, 3);
    L2->addWidget(new QLabel("Num points:", this), 0, 4);
    L2->addWidget(NPE = new QLineEdit("10000", this), 0, 5);
    L->addLayout(L2, 8, 0, 1, 8);
    L->setRowMinimumHeight(9, 20);
    QGridLayout *L3 = new QGridLayout();
    L3->addWidget(OK = new QPushButton("OK", this), 0, 0);
    L3->addWidget(Cancel = new QPushButton("Cancel", this), 0, 2);
    L->addLayout(L3, 10, 0, 1, 8);
    connect(addB, SIGNAL(clicked()), this, SLOT(addC()));
    connect(delB, SIGNAL(clicked()), this, SLOT(delC()));
    connect(setB, SIGNAL(clicked()), this, SLOT(setV()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

void CMLRPotDialog::addC()
{
    int n = LRCB->currentIndex();
    if ((*CVL)[n].isEmpty()) LRC->addItem(LRCB->currentText());
    else LRC->addItem(LRCB->currentText() + " = " + (*CVL)[n]);
}

void CMLRPotDialog::delC()
{
    if (LRC->currentRow() >= 0) delete LRC->takeItem(LRC->currentRow());
}

void CMLRPotDialog::setV()
{
    if (LRC->currentRow() < 0) return;
    QString T = LRC->currentItem()->text();
    int n = T.indexOf(" = ");
    if (n >= 0) T = T.left(n);
    if (CValE->text().isEmpty()) T += " = " + CValE->text();
    LRC->currentItem()->setText(T);
}

#endif
