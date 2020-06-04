//
// C++ Implementation: SFQSCalcDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "SFQSCalcDialog.h"


SFQSCalcDialog::SFQSCalcDialog(MainWindow* mw, QString Dir) : QDialog(mw)
{
    int n, N = (MW = mw)->getNumMolecules();
    setWindowTitle("Settings for calculation of FQS series:");
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Molecule:", this), 0, 0);
    L->addWidget(MolBox = new QComboBox(this), 0, 1, 1, 2);
    for (n=0; n<N; n++) MolBox->addItem(MW->getMolecule(n)->getName());
    MolBox->setEditable(false);
    L->addWidget(new QLabel("El. State:", this), 1, 0);
    L->addWidget(StateBox = new QComboBox(this), 1, 1, 1, 2);
    StateBox->setEditable(false);
    L->addWidget(new QLabel("Origin pot. of series:", this), 2, 0);
    L->addWidget(PotBox = new QComboBox(this), 2, 1, 1, 2);
    PotBox->setEditable(false);
    L->addWidget(new QLabel("Dir of potentials:", this), 3, 0);
    L->addWidget(PotDirE = new QLineEdit(Dir, this), 3, 1);
    L->addWidget(sPotDirB = new QPushButton("...", this), 3, 2);
    L->addWidget(new QLabel("Dir of fit datasets:", this), 4, 0);
    L->addWidget(FDatDirE = new QLineEdit(Dir, this), 4, 1);
    L->addWidget(sFDatDirB = new QPushButton("...", this), 4, 2);
    L->addWidget(new QLabel("Out.Turn.Point f. SFQS:", this), 5, 0);
    L->addWidget(SFQSRadE = new QLineEdit("12.0", this), 5, 1, 1, 2);
    L->addWidget(new QLabel("Num. paralell calc.:", this), 6, 0);
    L->addWidget(NumParFitsE = new QLineEdit("8", this), 6, 1, 1, 2);
    SFQSRadE->setValidator(new QDoubleValidator(SFQSRadE));
    NumParFitsE->setValidator(new QIntValidator(1, 100, NumParFitsE));
    L->setRowMinimumHeight(7, 20);
    L->addWidget(OK = new QPushButton("OK", this), 8, 0);
    L->addWidget(Cancel = new QPushButton("Cancel", this), 8, 1, 1, 2);
    connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(MolChanged(int)));
    connect(StateBox, SIGNAL(currentIndexChanged(int)), this, SLOT(StateChanged(int)));
    connect(sPotDirB, SIGNAL(clicked()), this, SLOT(sPotDir()));
    connect(sFDatDirB, SIGNAL(clicked()), this, SLOT(sFDatDir()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    MolChanged(0);
}

void SFQSCalcDialog::exec(ElState*& St, Potential*& OPot, QString& PotDir, QString& FDatDir, double& SFQSRad, int& NumParFits)
{
    if (QDialog::exec() == QDialog::Accepted)
    {
        St = MW->getMolecule(MolBox->currentIndex())->getStateP(StateBox->currentIndex());
        OPot = St->getPotential(PotBox->currentIndex());
        PotDir = PotDirE->text();
        FDatDir = FDatDirE->text();
        SFQSRad = SFQSRadE->text().toDouble();
        NumParFits = NumParFitsE->text().toInt();
    }
    else
    {
        St = 0;
        SFQSRad = 0.0;
        NumParFits = 0;
    }
}

void SFQSCalcDialog::MolChanged(int n)
{
    Molecule *Mol = MW->getMolecule(n);
    int N = Mol->getNumStates();
    StateBox->clear();
    for (n=0; n<N; n++) StateBox->addItem(Mol->getState(n));
    StateChanged(0);
}

void SFQSCalcDialog::sFDatDir()
{
    QString Res = QFileDialog::getExistingDirectory(this, "MolSpektAnalysis", FDatDirE->text());
    if (!Res.isEmpty()) FDatDirE->setText(Res);
}

void SFQSCalcDialog::sPotDir()
{
    QString Res = QFileDialog::getExistingDirectory(this, "MolSpektAnalysis", PotDirE->text());
    if (!Res.isEmpty()) PotDirE->setText(Res);
}

void SFQSCalcDialog::StateChanged(int n)
{
    ElState *St = (n>=0 ? MW->getMolecule(MolBox->currentIndex())->getStateP(n) : 0);
    PotBox->clear();
    if (St != 0) for (n=0; n < St->getNumPotentials(); n++) PotBox->addItem(St->getPotentialName(n));
}
