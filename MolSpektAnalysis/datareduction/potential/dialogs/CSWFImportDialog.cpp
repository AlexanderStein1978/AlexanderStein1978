//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "CSWFImportDialog.h"
#include "MainWindow.h"
#include "molecule.h"
#include "isotab.h"

#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QListWidget>
#include <QDir>
#include <QMessageBox>
#include <QFileDialog>


CSWFImportDialog::CSWFImportDialog(MainWindow* MW, QString InitDir): QDialog(MW)
{
    int n, N = MW->getNumMolecules();
    mw = MW;
    ST = CT = 0;
    setWindowTitle("MolSpektAnalysis");
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Import wave functions:", this), 0, 0, 1, 5);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Molecule:", this), 2, 0);
    L->addWidget(MolBox = new QComboBox(this), 2, 1);
    L->addWidget(new QLabel("Potential file:", this), 4, 0);
    L->addWidget(PotFile = new QLineEdit(InitDir, this), 4, 1);
    L->addWidget(SelectPotFile = new QPushButton("...", this), 4, 2);
    L->addWidget(new QLabel("File T5_E_Bv:", this), 5, 0);
    L->addWidget(T5_E_Bv = new QLineEdit(InitDir, this), 5, 1);
    L->addWidget(SelectT5_E_Bv = new QPushButton("...", this), 5, 2);
    L->addWidget(new QLabel("Electronic states:", this), 2, 3, 1, 2);
    L->addWidget(SList = new QListWidget(this), 3, 3, 3, 1);
    L->addWidget(SBox = new QComboBox(this), 3, 4);
    L->addWidget(AddState = new QPushButton("<- Add", this), 4, 4);
    L->addWidget(RemoveState = new QPushButton("<- Remove", this), 5, 4);
    L->setRowMinimumHeight(6, 20);
    L->addWidget(new QLabel("Data directories:", this), 7, 0, 1, 3);
    L->addWidget(new QLabel("Isotopologue:", this), 7, 3);
    L->addWidget(IsoBox = new QComboBox(this), 7, 4);
    L->addWidget(new QLabel("Directory for isotopologue:", this), 8, 0, 1, 3);
    L->addWidget(IsoDir = new QLineEdit(InitDir, this), 8, 3);
    L->addWidget(SelectIsoDir = new QPushButton("...", this), 8, 4);
    L->addWidget(IsoList = new QListWidget(this), 9, 0, 3, 4);
    L->addWidget(AddIso = new QPushButton("<- Add", this), 9, 4);
    L->setRowMinimumHeight(10, 20);
    L->addWidget(RemoveIso = new QPushButton("<- Remove", this), 11, 4);
    L->setRowMinimumHeight(12, 20);
    L->addWidget(OK = new QPushButton("OK", this), 13, 0, 1, 2);
    L->addWidget(Cancel = new QPushButton("Cancel", this), 13, 3, 1, 2);
    MolBox->setEditable(false);
    IsoBox->setEditable(false);
    SBox->setEditable(false);
    for (n=0; n<N; n++) MolBox->addItem(MW->getMolecule(n)->getName());
    connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(molChanged(int)));
    connect(SelectPotFile, SIGNAL(clicked()), this, SLOT(selectPotFile()));
    connect(SelectT5_E_Bv, SIGNAL(clicked()), this, SLOT(selectT5_E_Bv()));
    connect(SelectIsoDir, SIGNAL(clicked()), this, SLOT(selectIsoDir()));
    connect(AddIso, SIGNAL(clicked()), this, SLOT(addIso()));
    connect(AddState, SIGNAL(clicked()), this, SLOT(addState()));
    connect(RemoveIso, SIGNAL(clicked()), this, SLOT(removeIso()));
    connect(RemoveState, SIGNAL(clicked()), this, SLOT(removeState()));
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
    molChanged(0);
}

void CSWFImportDialog::addIso()
{
    QString Iso = IsoBox->currentText(), File = IsoDir->text();
    QString LI = (IsoList->count() > 0 ? IsoList->item(IsoList->count() - 1)->text() : "");
    QDir Dir(File);
    int n;
    if (!Iso.isEmpty() && ((n = LI.indexOf(' ')) > 0 ? File != LI.right(LI.length() - n - 1) : true) && Dir.exists())
    {
        n = IsoBox->currentIndex();
        IsoList->addItem(Iso + ' ' + File);
        IsoBox->removeItem(n);
        IL2.append(IL1.takeAt(n));
    }
    else QMessageBox::information(this, "MolSpektAnalysis", 
        "Error: You can add an isotopologue only if a valid isotopologue is choosen and a valid directory, which has to be different compared to the ones already choosen for other isotopologues!");
}

void CSWFImportDialog::addState()
{
    int n;
    if (!SBox->currentText().isEmpty())
    {
        SList->addItem(SBox->currentText());
        SBox->removeItem(n = SBox->currentIndex());
        SL2.append(SL1.takeAt(n));
    }
}

void CSWFImportDialog::getData(Molecule*& Mol, ElState**& States, int *&Component, int& NStates, int*& Iso, QString*& IsoDirs, 
                               int& NIso, QString &PotentialFile, QString &T5_E_BvFile)
{
    int n, m;
    QString IsoText;
    Mol = mol;
    States = new ElState*[NStates = SL2.count()];
    Component = new int[NStates];
    for (n=0; n < NStates; n++) 
    {
        States[n] = Mol->getStateP(ST[SL2[n]]);
        Component[n] = CT[SL2[n]];
    }
    Iso = new int[NIso = IL2.count()];
    IsoDirs = new QString[NIso];
    for (n=0; n < NIso; n++)
    {
        Iso[n] = IL2[n];
        IsoText = IsoList->item(n)->text();
        m = IsoText.indexOf(' ');
        IsoDirs[n] = IsoText.right(IsoText.length() - m - 1);
    }
    PotentialFile = PotFile->text();
    T5_E_BvFile = T5_E_Bv->text();
}

void CSWFImportDialog::molChanged(int i)
{
    int n, m, N, NC, M;
    IsoBox->clear();
    IsoList->clear();
    IL1.clear();
    IL2.clear();
    SL1.clear();
    SL2.clear();
    SList->clear();
    SBox->clear();
    mol = mw->getMolecule(i);
    if (mol == 0) return;
    IsoTab *Iso = mol->getIso();
    for (n = NC = 0, N = mol->getNumStates(); n < N; n++) 
        NC += (int(2.0 * mol->getStateP(n)->getS()) + 1) * (mol->getStateP(n)->getLambda() == 0 ? 1 : 2);
    if (ST != 0)
    {
        delete[] ST;
        delete[] CT;
    }
    ST = new int[NC];
    CT = new int[NC];
    for (n = NC = 0; n<N; n++)
    {
        M = int(2.0 * mol->getStateP(n)->getS()) + 1;
        if (M==1)
        {
            if (mol->getStateP(n)->getLambda() == 0)
            {
                ST[NC] = n;
                CT[NC] = 0;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n));
            }
            else
            {
                ST[NC] = n;
                CT[NC] = 0;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n) + " e");
                ST[NC] = n;
                CT[NC] = 1;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n) + " f");
            }
        }
        else for (m=0; m<M; m++)
        {
            if (mol->getStateP(n)->getLambda() == 0)
            {
                ST[NC] = n;
                CT[NC] = m;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n) + ' ' + QString::number(m));
            }
            else
            {
                ST[NC] = n;
                CT[NC] = 2*m;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n) + ' ' + QString::number(m) + 'e');
                ST[NC] = n;
                CT[NC] = 2*m+1;
                SL1.append(NC++);
                SBox->addItem(mol->getState(n) + ' ' + QString::number(m) + 'f');
            }
        }
    }
    for (n=0; n < Iso->numIso; n++)
    {
        IL1.append(n);
        IsoBox->addItem(Iso->getIsoName(n));
    }
    delete Iso;
}

void CSWFImportDialog::removeIso()
{
    QString Text = IsoList->currentItem()->text();
    int n = IsoList->currentRow();
    if (Text.isEmpty()) return;
    IL1.append(IL2.takeAt(n));
    delete IsoList->takeItem(n);
    IsoBox->addItem(Text.left(Text.indexOf(' ')));
}

void CSWFImportDialog::removeState()
{
    int n = SList->currentRow();
    if (n < 0) return;
    QListWidgetItem *Item = SList->takeItem(n);
    SBox->addItem(Item->text());
    delete Item;
    SL1.append(SL2.takeAt(n));
}

void CSWFImportDialog::selectIsoDir()
{
    QString Dir = QFileDialog::getExistingDirectory(this, "Select directory with wave functions for the isotopologue " 
                                                    + IsoBox->currentText(), IsoDir->text());
    if (!Dir.isEmpty()) IsoDir->setText(Dir);
}

void CSWFImportDialog::selectPotFile()
{
    QString File = QFileDialog::getOpenFileName(this, "Please select the potential file", PotFile->text(), "Potential files (*)");
    if (!File.isEmpty()) PotFile->setText(File);
}

void CSWFImportDialog::selectT5_E_Bv()
{
    QString File = QFileDialog::getOpenFileName(this, "Please select the file T5_E_Bv", T5_E_Bv->text(), "T5_E_Bv (T5_E_Bv)");
    if (!File.isEmpty()) T5_E_Bv->setText(File);
}
