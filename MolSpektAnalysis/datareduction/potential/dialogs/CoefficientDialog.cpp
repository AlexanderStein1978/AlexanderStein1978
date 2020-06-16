//
// C++ Implementation: CoefficientDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "CoefficientDialog.h"

#include <QListWidget>
#include <QLineEdit>
#include <QGridLayout>
#include <QPushButton>
#include <QLabel>


CoefficientDialog::CoefficientDialog(QWidget *parent, QStringList &Data) : QDialog(parent)
{
    QGridLayout *Layout = new QGridLayout(this);
    Add = new QPushButton("Add", this);
    Layout->addWidget(Add, 1, 1);
    connect(Add, SIGNAL(clicked()), this, SLOT(add()));
    Cancel = new QPushButton("Cancel", this);
    Layout->addWidget(Cancel, 4, 1);
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    CList = new QListWidget(this);
    CList->addItems(Data);
    Layout->addWidget(CList, 1, 0, 2, 1);
    Edit = new QLineEdit(this);
    Layout->addWidget(Edit, 0, 1);
    Label = new QLabel("Name/R [A]:", this);
    Layout->addWidget(Label, 0, 0);
    OK = new QPushButton("OK", this);
    Layout->addWidget(OK, 4, 0);
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    Remove = new QPushButton("Remove", this);
    Layout->addWidget(Remove, 2, 1);
    connect(Remove, SIGNAL(clicked()), this, SLOT(Remove()));
    Layout->setRowMinimumHeight(3, 10);
}

void CoefficientDialog::add()
{
    QString C = Edit->text();
    if (!C.isEmpty()) CList->addItem(C);
}

QStringList CoefficientDialog::getResults()
{
    int n, N = CList->count();
    QStringList Result;
    for (n=0; n<N; n++) Result << CList->item(n)->text();
    return Result;
}

void CoefficientDialog::remove()
{
    int n = CList->currentRow();
    if (n<0) return;
    QListWidgetItem *CI = CList->takeItem(n);
    delete CI;
}
