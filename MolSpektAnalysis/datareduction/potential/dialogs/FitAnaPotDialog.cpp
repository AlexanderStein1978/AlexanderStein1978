//
// C++ Implementation: FitAnaPotDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "FitAnaPotDialog.h"

#include <QPushButton>
#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QListWidget>


FitAnaPotDialog::FitAnaPotDialog(bool improve, double Ri, double Ra, double iExp, int NC, int NLRC, int* pLRC, double hCi_cP, 
                                 PotentialType CurrentPotential, QWidget* parent) : QDialog(parent)
{
    int n;
    QPushButton *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    QGridLayout *L = new QGridLayout(this);
    setWindowTitle("Fit potential");
    L->addWidget(new QLabel("Fit potential of type:", this), 0, 0);
    L->addWidget(TypeBox = new QComboBox(this), 0, 1);
    TypeBox->addItems(QStringList() << "spline potential" << "analytical potential");
    TypeBox->setEditable(false);
    switch (CurPotential = CurrentPotential)
    {
        case SplinePotential:
            TypeBox->setCurrentIndex(0);
            break;
        case analyticalPotential:
        default:
            TypeBox->setCurrentIndex(1);
            break;
    }
    L->addWidget(Improve = new QCheckBox("Fit current potential", this), 1, 0, 1, 2);
    Improve->setChecked(improve);
    Improve->setEnabled(improve);
    L->addWidget(UseSVD = new QCheckBox("Use Singular Value Decomposition", this), 2, 0, 1, 2);
    UseSVD->setChecked(true);
    L->addWidget(UseLeveMrq = new QCheckBox("Use Levenberg-Marquardt", this), 3, 0, 1, 2);
    UseLeveMrq->setChecked(true);
    L->addWidget(AdjustTEWF = new QCheckBox("Adjust TE weight factors", this), 4, 0, 1, 2);
    AdjustTEWF->setChecked(false);
    L->addWidget(new QLabel("Ri:", this), 5, 0);
    L->addWidget(RiE = new QLineEdit(QString::number(Ri, 'f', 3), this), 5, 1);
    L->addWidget(new QLabel("Ra:", this), 6, 0);
    L->addWidget(RaE = new QLineEdit(QString::number(Ra, 'f', 1), this), 6, 1);
    L->addWidget(new QLabel("iExp:", this), 7, 0);
    L->addWidget(iExpE = new QLineEdit(QString::number(iExp, 'g', 12), this), 7, 1);
    iExpE->setValidator(new QDoubleValidator(2.0, 20.0, 12, iExpE));
    L->addWidget(new QLabel("Number of coefficients:", this), 8, 0);
    L->addWidget(NCE = new QLineEdit(QString::number(NC >= 4 ? NC : 15), this), 8, 1);
    NCE->setValidator(new QIntValidator(4, 40, NCE));
    L->setRowMinimumHeight(9, 20);
    L->addWidget(new QLabel("Long range coefficients:", this), 10, 0, 1, 2);
    L->addWidget(CoeffList = new QListWidget(this), 11, 0, 3, 1);
    for (n=0; n < NLRC; n++) CoeffList->addItem('C' + QString::number(pLRC[n]));
    L->addWidget(CoeffBox = new QComboBox(this), 11, 1);
    CoeffBox->addItems(QStringList() << "C3" << "C6" << "C8" << "C10" << "C12");
    CoeffBox->setEditable(true);
    L->addWidget(Add = new QPushButton("<- Add", this), 12, 1);
    L->addWidget(Remove = new QPushButton("Remove", this), 13, 1);
    L->addWidget(new QLabel("high Ci rel. prec.:", this), 14, 0);
    L->addWidget(hCi_cPE = new QLineEdit(QString::number(hCi_cP, 'g', 4), this), 14, 1);
    if (hCi_cP < 0.0) hCi_cPE->setEnabled(false);
    L->setRowMinimumHeight(15, 20);
    L->addWidget(new QLabel("Fit range:", this), 16, 0);
    L->addWidget(new QLabel("R min:", this), 17, 0);
    L->addWidget(RminE = new QLineEdit(QString::number(Ri - 1.0, 'f', 3), this), 17, 1);
    L->addWidget(new QLabel("R max:", this), 18, 0);
    L->addWidget(RmaxE = new QLineEdit("50.0", this), 18, 1);
    L->addWidget(new QLabel("Number of points:", this), 19, 0);
    L->addWidget(NumPE = new QLineEdit("10000", this), 19, 1);
    L->addWidget(new QLabel("b min:", this), 20, 0);
    L->addWidget(bminE = new QLineEdit("-0.1", this), 20, 1);
    L->addWidget(new QLabel("b max:", this), 21, 0);
    L->addWidget(bmaxE = new QLineEdit("-0.59", this), 21, 1);
    L->addWidget(new QLabel("b step:", this), 22, 0);
    L->addWidget(bstepE = new QLineEdit("-0.01", this), 22, 1);
    L->setRowMinimumHeight(23, 20);
    L->addWidget(OK, 24, 0);
    L->addWidget(Cancel, 24, 1);
    connect(TypeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(typeChanged(int)));
    connect(Improve, SIGNAL(toggled(bool)), this, SLOT(improveChanged(bool)));
    connect(Add, SIGNAL(clicked()), this, SLOT(addCoeff()));
    connect(Remove, SIGNAL(clicked()), this, SLOT(removeCoeff()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    improveChanged(improve);
    typeChanged(TypeBox->currentIndex());
}

void FitAnaPotDialog::addCoeff()
{
    QString C = CoeffBox->currentText();
    if (C[0] == 'c' || C[0] == 'C') C = C.right(C.length() - 1);
    int n, c = C.toInt();
    if (c>0)
    {
        C = 'C' + QString::number(c);
        for (n=0; (n < CoeffList->count() ? CoeffList->item(n)->text() != C : false); n++) ;
        if (n == CoeffList->count()) CoeffList->addItem(C);
    }
}

bool FitAnaPotDialog::exec(bool &improve, bool &useLeveMrq, bool &useSVD, double& Ri, double& Ra, double& iExp, int& NC, int& NLRC, 
                           int*& pLRC, double& Rmin, double& Rmax, int& NumP, double &bmin, double &bmax, double &bstep, bool &ATEWF,
                           double &hCi_cP, PotentialType &NewPotType)
{
    if (QDialog::exec() == QDialog::Rejected) return false;
    improve = Improve->isChecked();
    useLeveMrq = UseLeveMrq->isChecked();
    useSVD = UseSVD->isChecked();
    ATEWF = AdjustTEWF->isChecked();
    hCi_cP = hCi_cPE->text().toDouble();
    if (improve) return true;
    int n;
    QString C;
    Ri = RiE->text().toDouble();
    Ra = RaE->text().toDouble();
    iExp = iExpE->text().toDouble();
    NC = NCE->text().toDouble();
    pLRC = new int[NLRC = CoeffList->count()];
    for (n=0; n < NLRC; n++) 
    {
        C = CoeffList->item(n)->text();
        pLRC[n] = C.right(C.length() - 1).toInt();
    }
    Rmin = RminE->text().toDouble();
    Rmax = RmaxE->text().toDouble();
    NumP = NumPE->text().toInt();
    bmin = bminE->text().toDouble();
    bmax = bmaxE->text().toDouble();
    bstep = bstepE->text().toDouble();
    NewPotType = (TypeBox->currentIndex() == 0 ? SplinePotential : analyticalPotential);
    return true;
}

void FitAnaPotDialog::improveChanged(bool State)
{
    Add->setEnabled(!State);
    CoeffBox->setEnabled(!State);
    CoeffList->setEnabled(!State);
    NCE->setEnabled(!State);
    NumPE->setEnabled(!State);
    RaE->setEnabled(!State);
    Remove->setEnabled(!State);
    RiE->setEnabled(!State);
    RmaxE->setEnabled(!State);
    RminE->setEnabled(!State);
    iExpE->setEnabled(!State);
    bminE->setEnabled(!State);
    bmaxE->setEnabled(!State);
    bstepE->setEnabled(!State);
}

void FitAnaPotDialog::removeCoeff()
{
    QListWidgetItem *I = CoeffList->takeItem(CoeffList->currentRow());
    delete I;
}

void FitAnaPotDialog::typeChanged(int Index)
{
    if ((CurPotential == SplinePotential && Index == 0) || (CurPotential == analyticalPotential && Index == 1)) 
        Improve->setEnabled(true);
    else
    {
        if (Improve->isChecked()) Improve->setChecked(false);
        Improve->setEnabled(false);
    }
}
