//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.h"
#include "elstate.h"
#include "MainWindow.h"
#include "molecule.h"

#include <QGridLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>
#include <QFileDialog>


CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::CreateAnaPotSeriesFromMCSSplinePotSeriesDialog(MainWindow* parent)
                                                : QDialog(parent), MW(parent)
{
    int n;
    QGridLayout *L = new QGridLayout(this);
    QPushButton *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    QPushButton *FDDB = new QPushButton("...", this), *APDB = new QPushButton("...", this);
    setWindowTitle("Create AnaPot series from MCS SplinePot series");
    L->addWidget(improve = new QCheckBox("Improve existing analytical potentials", this), 0, 0, 1, 3);
    L->addWidget(new QLabel("Dir with SplinePots:", this), 1, 0);
    L->addWidget(SplinePotDirE = new QLineEdit(MW->getDir(MDIChild::PotData), this), 1, 1);
    L->addWidget(SPDB = new QPushButton("...", this), 1, 2);
    L->addWidget(new QLabel("Dir with FitData sets:", this), 2, 0);
    L->addWidget(FDDirE = new QLineEdit(MW->getDir(MDIChild::FitDataSet), this), 2, 1);
    L->addWidget(FDDB, 2, 2);
    L->addWidget(new QLabel("Dir for AnaPots:", this), 3, 0);
    L->addWidget(AnaPotDirE = new QLineEdit(MW->getDir(MDIChild::PotData), this), 3, 1);
    L->addWidget(APDB, 3, 2);
    L->addWidget(new QLabel("Molecule:", this), 4, 0);
    L->addWidget(MoleculeBox = new QComboBox(this), 4, 1, 1, 2);
    L->addWidget(new QLabel("Electronic state:", this), 5, 0);
    L->addWidget(StateBox = new QComboBox(this), 5, 1, 1, 2);
    L->addWidget(new QLabel("Fit precision:", this), 6, 0);
    L->addWidget(PrecE = new QLineEdit("1.0", this), 6, 1, 1, 2);
    L->addWidget(new QLabel("Max N fit iterations:", this), 7, 0);
    L->addWidget(MaxIt = new QComboBox(this), 7, 1, 1, 2);
    L->addWidget(useSVD = new QCheckBox("Use SVD", this), 8, 0);
    L->addWidget(useLM = new QCheckBox("Use Levenb.Mrq.", this), 8, 1, 1, 2);
    L->setRowMinimumHeight(9, 20);
    L->addWidget(OK, 10, 0);
    L->addWidget(Cancel, 10, 1, 1, 2);
    for (n=0; n < MW->getNumMolecules(); n++) MoleculeBox->addItem(MW->getMolecule(n)->getName());
    MoleculeBox->setEditable(false);
    StateBox->setEditable(false);
    PrecE->setValidator(new QDoubleValidator(1e-10, 1e10, 10, PrecE));
    improve->setChecked(false);
    MaxIt->addItem("No limit");
    for (n=1; n<=10; n++) MaxIt->addItem(QString::number(n));
    useLM->setChecked(true);
    useSVD->setChecked(true);
    connect(improve, SIGNAL(toggled(bool)), this, SLOT(improveChanged(bool)));
    connect(MaxIt, SIGNAL(editTextChanged(QString)), this, SLOT(MaxItChanged(QString)));
    connect(SPDB, SIGNAL(clicked()), this, SLOT(showSplinePotDirDialog()));
    connect(FDDB, SIGNAL(clicked()), this, SLOT(showFDDirDialog()));
    connect(APDB, SIGNAL(clicked()), this, SLOT(showAnaPotDirDialog()));
    connect(MoleculeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(molChanged(int)));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    molChanged(0);
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::getData(QString& MolFN, QString& StateN, QString& FDDir,
            QString& SPDir, QString& APDir, bool& improveAnaPots, bool& UseSvd, bool& UseLeveMarq, int& MaxIter, 
            double& Prec)
{
    MolFN = MW->getMolecule(MoleculeBox->currentIndex())->getFileName();
    StateN = StateBox->currentText();
    FDDir = FDDirE->text();
    SPDir = SplinePotDirE->text();
    APDir = AnaPotDirE->text();
    improveAnaPots = improve->isChecked();
    UseSvd = useSVD->isChecked();
    UseLeveMarq = useLM->isChecked();
    MaxIter = MaxIt->currentText().toInt();
    Prec = PrecE->text().toDouble();
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::improveChanged(bool checked)
{
    SplinePotDirE->setEnabled(!checked);
    SPDB->setEnabled(!checked);
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::MaxItChanged(const QString& text)
{
    if (text.toInt() <= 0) MaxIt->setCurrentIndex(0);
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::molChanged(int i)
{
    int n;
    Molecule *mol = MW->getMolecule(i);
    StateBox->clear();
    for (n=0; n < mol->getNumStates(); n++) StateBox->addItem(mol->getState(n));
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::showAnaPotDirDialog()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select directory for AnaPots", AnaPotDirE->text());
    if (!Dir.isEmpty()) AnaPotDirE->setText(Dir);
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::showFDDirDialog()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select directory with FitData sets", FDDirE->text());
    if (!Dir.isEmpty()) FDDirE->setText(Dir);
}

void CreateAnaPotSeriesFromMCSSplinePotSeriesDialog::showSplinePotDirDialog()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select directory with SplinePots", SplinePotDirE->text());
    if (!Dir.isEmpty()) SplinePotDirE->setText(Dir);
}
