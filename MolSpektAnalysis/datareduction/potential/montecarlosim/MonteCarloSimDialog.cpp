//
// C++ Implementation: MonteCarloSimDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "MonteCarloSimDialog.h"
#include "MainWindow.h"

#include <QGridLayout>
#include <QPushButton>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QLabel>
#include <QFileInfo>
#include <QFileDialog>
#include <QMessageBox>


MonteCarloSimDialog::MonteCarloSimDialog(MainWindow* MW) : QDialog(MW)
{
    QGridLayout *L = new QGridLayout(this);
    QPushButton *FDB = new QPushButton("...", this), *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
    QDoubleValidator *Val;
    setWindowTitle("Monte Carlo Simulation");
    L->addWidget(new QLabel("Number of iterations:", this), 0, 0);
    L->addWidget(NItE = new QLineEdit("1000", this), 0, 1, 1, 2);
    L->addWidget(new QLabel("Number of parallel fits:", this), 1, 0);
    L->addWidget(NParFitsE = new QLineEdit("8", this), 1, 1, 1, 2);
    L->addWidget(new QLabel("Result directory:", this), 2, 0);
    L->addWidget(DirectoryE = new QLineEdit(MW->getDir(), this), 2, 1);
    L->addWidget(FDB, 2, 2);
    L->addWidget(new QLabel("Uncertainty factor:", this), 3, 0);
    L->addWidget(UncFactE = new QLineEdit("1.0", this), 3, 1, 1, 2);
    L->setRowMinimumHeight(4, 20);
    L->addWidget(OK, 5, 0);
    L->addWidget(Cancel, 5, 1, 1, 2);
    NItE->setValidator(new QIntValidator(1, 10000, NItE));
    NParFitsE->setValidator(new QIntValidator(1, 100, NParFitsE));
    UncFactE->setValidator(Val = new QDoubleValidator(UncFactE));
    Val->setBottom(1e-300);
    connect(FDB, SIGNAL(clicked()), this, SLOT(showFileDialog()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

bool MonteCarloSimDialog::exec(QString& Directory, int& NIterations, int& NParFit, double &UncFact)
{
    if (QDialog::exec() == QDialog::Rejected) return false;
    QFileInfo FI(Directory = DirectoryE->text());
    if (!FI.isDir())
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Error: " + Directory + " is not an existing directory!");
        return exec(Directory, NIterations, NParFit, UncFact);
    }
    if (Directory.right(1) != "\\" && Directory.right(1) != "/") Directory += DIRSEP;
    NParFit = NParFitsE->text().toInt();
    NIterations = NItE->text().toInt();
    UncFact = UncFactE->text().toDouble();
    return true;
}

void MonteCarloSimDialog::showFileDialog()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select directory for results", DirectoryE->text());
    if (!Dir.isEmpty()) DirectoryE->setText(Dir);
}
