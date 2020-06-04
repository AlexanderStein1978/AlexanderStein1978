//
// C++ Implementation: ImprovePotSeriesDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "ImprovePotSeriesDialog.h"


ImprovePotSeriesDialog::ImprovePotSeriesDialog(MainWindow* MW, QString Dir, Potential* OPot): QDialog()
{
    QGridLayout *L = new QGridLayout(this);
    QPushButton *sPotDirB = new QPushButton("...", this), *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    QDoubleValidator *ThV;
    Molecule *TM;
    int n, N = (mw = MW)->getNumMolecules();
    Mol = (OPot != 0 ? OPot->getMolecule() : 0);
    St = (OPot != 0 ? OPot->getElState() : 0);
    setWindowTitle("Settings for Improvement of Pot Series");
    L->addWidget(new QLabel("Molecule:", this), 0, 0);
    L->addWidget(MoleculeB = new QComboBox(this), 0, 1, 1, 2);
    for (n=0; n<N; n++) 
    {
        MoleculeB->addItem((TM = MW->getMolecule(n))->getName());
        if (Mol == TM) MoleculeB->setCurrentIndex(n);
    }
    MoleculeB->setEditable(false);
    L->addWidget(new QLabel("Electronic State:", this), 1, 0);
    L->addWidget(StateB = new QComboBox(this), 1, 1, 1, 2);
    StateB->setEditable(false);
    L->addWidget(new QLabel("Fit Data set:", this), 2, 0);
    L->addWidget(FitDataB = new QComboBox(this), 2, 1, 1, 2);
    FitDataB->setEditable(false);
    L->addWidget(new QLabel("Directory of Potentials:"), 3, 0);
    L->addWidget(PotDirE = new QLineEdit(Dir, this), 3, 1);
    L->addWidget(sPotDirB, 3, 2);
    L->addWidget(new QLabel("Threshol:", this), 4, 0);
    L->addWidget(ThresholdE = new QLineEdit((OPot != 0 ? (OPot->getPotType() == SplinePotential ? "0.0" : "1.0") : "1.0"), this), 
                 4, 1, 1, 2);
    ThresholdE->setValidator(ThV = new QDoubleValidator(ThresholdE));
    ThV->setBottom(0.0);
    L->addWidget(new QLabel("Number of parallel fits:", this), 5, 0);
    L->addWidget(NumParFitsE = new QLineEdit("8", this), 5, 1, 1, 2);
    NumParFitsE->setValidator(new QIntValidator(1, 100, NumParFitsE));
    L->setRowMinimumHeight(6, 20);
    L->addWidget(OK, 7, 0);
    L->addWidget(Cancel, 7, 1, 1, 2);
    connect(MoleculeB, SIGNAL(currentIndexChanged(int)), this, SLOT(MolChanged(int)));
    connect(StateB, SIGNAL(currentIndexChanged(int)), this, SLOT(StateChanged(int)));
    connect(sPotDirB, SIGNAL(clicked()), this, SLOT(sPotDir()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    MolChanged(0);
}

bool ImprovePotSeriesDialog::exec(ElState*& St, FitData*& FDat, QString& PotDir, double& Threshold, int &NumParFits)
{
    if (QDialog::exec() == QDialog::Rejected) return false;
    Molecule *Mol = mw->getMolecule(MoleculeB->currentIndex());
    St = Mol->getStateP(StateB->currentIndex());
    if (St == 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "You have to select a molecule where fit data is available!");
        return exec(St, FDat, PotDir, Threshold, NumParFits);
    }
    FDat = St->getFitData(FitDataB->currentIndex());
    if (FDat == 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "You have to select an electronic state where fit data is available!");
        return exec(St, FDat, PotDir, Threshold, NumParFits);
    }
    QDir D(PotDir = PotDirE->text());
    QStringList L = D.entryList(QStringList("*.pot"), QDir::Files);
    if (L.count() == 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "You have to select a containing potential files (\"*.pot\")!");
        return exec(St, FDat, PotDir, Threshold, NumParFits);
    }
    Threshold = ThresholdE->text().toDouble();
    NumParFits = NumParFitsE->text().toInt();
    return true;
}

void ImprovePotSeriesDialog::MolChanged(int i)
{
    Mol = mw->getMolecule(i);
    int n, N = Mol->getNumStates();
    StateB->clear();
    for (n=0; n<N; n++) 
    {
        StateB->addItem(Mol->getState(n));
        if (Mol->getStateP(n) == St) StateB->setCurrentIndex(n);
    }
}

void ImprovePotSeriesDialog::sPotDir()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select a directory containing potentials (\"*.pot\")", PotDirE->text());
    if (!Dir.isEmpty()) PotDirE->setText(Dir);
}

void ImprovePotSeriesDialog::StateChanged(int i)
{
    St = Mol->getStateP(i);
    FitDataB->clear();
    if (St != 0)
    {
        int n, N = St->getNumFitDataSets();
        for (n=0; n<N; n++) FitDataB->addItem(St->getFitDataName(n));
        FitDataB->setCurrentIndex(St->getMainFitDataNum());
    }
}
