//
// C++ Implementation: vAssignDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "vassigndialog.h"
#include "MainWindow.h"
#include "elstate.h"
#include "molecule.h"

#include <QGridLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>


vAssignDialog::vAssignDialog(MainWindow* parent, ElState* S) : QDialog(parent)
{
    int n, N = (State = S)->getNumTermTables();
    QGridLayout *L = new QGridLayout(this);
    setWindowTitle("Please select data source for v assignment");
    L->addWidget(new QLabel("Electronic state:", this), 0, 0);
    L->addWidget(StateBox = new QComboBox(this), 0, 1);
    L->addWidget(new QLabel("Data Source:", this), 1, 0);
    L->addWidget(Box = new QComboBox(this), 1, 1);
    L->addWidget(List = new QListWidget(this), 2, 0, 3, 1);
    QGridLayout *L2 = new QGridLayout();
    L2->addWidget(new QLabel("v=", this), 0, 0);
    L2->addWidget(vMin = new QLineEdit(this), 0, 1);
    L2->addWidget(new QLabel("to", this), 0, 2);
    L2->addWidget(vMax = new QLineEdit(this), 0, 3);
    L->addLayout(L2, 2, 1, 1, 1);
    L->addWidget(Add = new QPushButton("<-Add", this), 3, 1);
    L->setRowMinimumHeight(4, 20);
    L->addWidget(Remove = new QPushButton("<-Remove", this), 4, 1);
    L->addWidget(new QLabel("Tolerance:", this), 5, 0);
    L->addWidget(Tolerance = new QLineEdit("10.0", this), 5, 1);
    L->addWidget(new QLabel("Tol f. dbl. obs.:", this), 6, 0);
    L->addWidget(DATolerance = new QLineEdit("0.03", this), 6, 1);
    L->setRowMinimumHeight(7, 20);
    L->addWidget(OK = new QPushButton("OK", this), 8, 0);
    L->addWidget(Cancel = new QPushButton("Cancel", this), 8, 1);
    for (n=0; n < N; n++) Box->addItem(S->getTermTableName(n));
    Box->setEditable(false);
    Molecule *Mol = State->getMolecule();
    for (n=0, N = Mol->getNumStates(); n<N; ++n)
    {
        StateBox->addItem(Mol->getState(n));
        if (Mol->getStateP(n) == State) StateBox->setCurrentIndex(n);
    }
    StateBox->setEditable(false);
    connect(StateBox, SIGNAL(currentIndexChanged(int)), this, SLOT(stateBoxChanged(int)));
    connect(Add, SIGNAL(clicked()), this, SLOT(addSource()));
    connect(Remove, SIGNAL(clicked()), this, SLOT(removeSource()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

void vAssignDialog::addSource()
{
    int n = Box->currentIndex();
    QString Buffer = Box->currentText();
    if (!vMin->text().isEmpty() || !vMax->text().isEmpty())
    {
        int vm = vMin->text().toInt(), vM = vMax->text().toInt();
        if (vm < 0) vm = 0;
        if (vM < vm) vM = -1;
        Buffer += " v=" + QString::number(vm) + "to" + QString::number(vM);
    }
    List->addItem(Buffer);
    Box->removeItem(n);
}

bool vAssignDialog::getSelectedTables(TermTable**& T, int *&vMin, int *&vMax, int& N, double &Tol, double &DATol)
{
    int n, i, j;
    if ((N = List->count()) == 0) return false;
    T = new TermTable*[N];
    vMin = new int[N];
    vMax = new int[N];
    Tol = Tolerance->text().toDouble();
    DATol = DATolerance->text().toDouble();
    Molecule *Mol = State->getMolecule();
    for (n=0; n<N; n++)
    {
        QString Buffer = List->item(n)->text();
        if ((i = Buffer.indexOf("v=")) > 0)
        {
            j = Buffer.indexOf("to");
            vMin[n] = Buffer.mid(i+2, j-i-2).toInt();
            vMax[n] = Buffer.right(Buffer.length() - j - 2).toInt();
            Buffer = Buffer.left(i-1);
        }
        else vMin[n] = vMax[n] = -1;
        T[n] = Mol->getTermTable(Buffer);
    }
    return true;
}

void vAssignDialog::removeSource()
{
    int n = List->currentRow();
    if (n<0) return;
    Box->addItem(List->currentItem()->text());
    QListWidgetItem* Item = List->takeItem(n);
    QString Buffer = Item->text();
    n = Buffer.indexOf("v=");
    if (n>0) Buffer = Buffer.left(n-1);
    if (State->getMolecule()->getStateOfTermTable(Buffer) == State) Box->addItem(Buffer);
    delete Item;
}

void vAssignDialog::stateBoxChanged(int index)
{
    State = State->getMolecule()->getStateP(index);
    Box->clear();
    int n, m, L = List->count(), N = State->getNumTermTables();
    QString TName;
    for (n=0; n<N; ++n)
    {
        for (m=0, TName = State->getTermTableName(n); m<L && List->item(m)->text() != TName; ++m) ;
        if (m==L) Box->addItem(TName);
    }
}
