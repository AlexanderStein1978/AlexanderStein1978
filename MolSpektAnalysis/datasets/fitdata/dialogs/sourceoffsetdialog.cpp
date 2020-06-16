//
// C++ Implementation: SourceOffsetDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "sourceoffsetdialog.h"
#include "MainWindow.h"

#include <QGridLayout>
#include <QComboBox>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>


SourceOffsetDialog::SourceOffsetDialog(QStringList& Names, double* Offsets, MainWindow* parent): QDialog(parent)
{
    int n;
    offsets = Offsets;
    lastIndex = 0;
    QGridLayout *L = new QGridLayout(this);
    QComboBox *NameB = new QComboBox(this);
    QPushButton *SaveB = new QPushButton("Save", this), *Abort = new QPushButton("Abort", this);
    setWindowTitle("Change source offset");
    L->addWidget(new QLabel("Source:", this), 0, 0);
    L->addWidget(NameB, 0, 1);
    NameB->setEditable(false);
    for (n=0; n < Names.count(); n++) NameB->addItem(Names[n]);
    L->addWidget(new QLabel("Offset:", this), 1, 0);
    L->addWidget(OffsetE = new QLineEdit(QString::number(Offsets[0], 'f', 12), this), 1, 1);
    L->setRowMinimumHeight(2, 20);
    L->addWidget(SaveB, 3, 0);
    L->addWidget(Abort, 3, 1);
    connect(NameB, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangeIndex(int)));
    connect(SaveB, SIGNAL(clicked()), this, SLOT(Save()));
    connect(Abort, SIGNAL(clicked()), this, SLOT(reject()));
}

void SourceOffsetDialog::ChangeIndex(int Index)
{
    offsets[lastIndex] = OffsetE->text().toDouble();
    OffsetE->setText(QString::number(offsets[lastIndex = Index], 'f', 12));
}

void SourceOffsetDialog::Save()
{
    offsets[lastIndex] = OffsetE->text().toDouble();
    accept();
}
