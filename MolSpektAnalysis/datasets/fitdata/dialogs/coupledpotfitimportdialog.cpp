//
// C++ Implementation: CoupledPotFitOutputImportDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "coupledpotfitimportdialog.h"


CoupledPotFitOutputImportDialog::CoupledPotFitOutputImportDialog(MainWindow *mw)
    : QDialog(mw)
{
    OKButton = new QPushButton("OK", this);
    CancelButton = new QPushButton("Cancel", this);
    AddButton = new QPushButton("Add", this);
    RemoveButton = new QPushButton("Remove", this);
    ChoosenFitDataSets = new QListWidget(this);
    MolBox = new QComboBox(this);
    StateBox = new QComboBox(this);
    FitDataBox = new QComboBox(this);
    Update = new QCheckBox("Update deviations in existing fitdata", this);
    MW = mw;
    Mol = 0;
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Molecule:", this), 0, 0);
    L->addWidget(MolBox, 0, 1);
    L->addWidget(Update, 1, 0, 1, 2);
    L->addWidget(new QLabel("v assign tolerancy [cm^-1]", this), 2, 0);
    L->addWidget(TolerancyEdit = new QLineEdit("0.2", this), 2, 1);
    L->setRowMinimumHeight(3, 20);
    L->addWidget(new QLabel("Electronic state:", this), 4, 0);
    L->addWidget(StateBox, 4, 1);
    L->addWidget(new QLabel("Choosen FitData sets:", this), 5, 0, 1, 2);
    L->addWidget(ChoosenFitDataSets, 6, 0, 3, 1);
    L->addWidget(FitDataBox, 6, 1);
    L->addWidget(AddButton, 7, 1);
    L->addWidget(RemoveButton, 8, 1);
    L->setRowMinimumHeight(9, 20);
    L->addWidget(OKButton, 10, 0);
    L->addWidget(CancelButton, 10, 1);
    Update->setChecked(true);
    TolerancyEdit->setValidator(new QDoubleValidator(0.0, 1E6, 50, TolerancyEdit));
    int n, N = MW->getNumMolecules();
    for (n=0; n<N; ++n) MolBox->addItem(MW->getMolecule(n)->getName());
    MolBox->setEditable(false);
    StateBox->setEditable(false);
    FitDataBox->setEditable(false);
    connect(OKButton, SIGNAL(clicked()), this, SLOT(OK()));
    connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));
    connect(AddButton, SIGNAL(clicked()), this, SLOT(Add()));
    connect(RemoveButton, SIGNAL(clicked()), this, SLOT(Remove()));
    connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(MoleculeChanged(int)));
    connect(StateBox, SIGNAL(currentIndexChanged(int)), this, SLOT(StateChanged(int)));
    connect(Update, SIGNAL(toggled(bool)), this, SLOT(UpdateChanged(bool)));
    MoleculeChanged(0);
}

void CoupledPotFitOutputImportDialog::getResult(Molecule *&mol, FitData **&fDats, int &NfDats, double &Tolerancy)
{
    int n, m = 0, N = ChoosenFitDataSets->count();
    mol = Mol;
    if (Mol != 0 && Update->isChecked())
    {
        fDats = new FitData*[N];
        for (n=m=0; n<N; ++n)
        {
            fDats[m] = Mol->getFitData(ChoosenFitDataSets->item(m)->text());
            if (fDats[m] == 0) delete ChoosenFitDataSets->takeItem(m);
            else ++m;
        }
    }
    else fDats = 0;
    NfDats = m;
    Tolerancy = TolerancyEdit->text().toDouble();
}

void CoupledPotFitOutputImportDialog::MoleculeChanged(int Index)
{
    Mol = MW->getMolecule(Index);
    StateBox->clear();
    if (Mol != 0)
    {
        int n, N = Mol->getNumStates();
        for (n=0; n<N; n++) StateBox->addItem(Mol->getState(n));
    }
}

void CoupledPotFitOutputImportDialog::UpdateChanged(bool update)
{
    TolerancyEdit->setEnabled(update);
    StateBox->setEnabled(update);
    ChoosenFitDataSets->setEnabled(update);
    FitDataBox->setEnabled(update);
    AddButton->setEnabled(update);
    RemoveButton->setEnabled(update);
}

void CoupledPotFitOutputImportDialog::StateChanged(int Index)
{
    FitDataBox->clear();
    if (Mol != 0)
    {
        ElState *S = Mol->getStateP(Index);
        if (S != 0)
        {
            int n, N = S->getNumFitDataSets();
            for (n=0; n<N; n++)
            {
                QString Name = S->getFitDataName(n);
                if (ChoosenFitDataSets->findItems(Name, Qt::MatchFixedString | Qt::MatchCaseSensitive).isEmpty()) FitDataBox->addItem(Name);
            }
        }
    }
}

void CoupledPotFitOutputImportDialog::Add()
{
    QString Name = FitDataBox->currentText();
    if (!Name.isEmpty())
    {
        ChoosenFitDataSets->addItem(Name);
        FitDataBox->removeItem(FitDataBox->currentIndex());
    }
}

void CoupledPotFitOutputImportDialog::Remove()
{
    QListWidgetItem* Item = ChoosenFitDataSets->takeItem(ChoosenFitDataSets->currentRow());
    if (Item != 0)
    {
        if (Mol != 0)
        {
            ElState *S = Mol->getStateP(StateBox->currentIndex());
            if (S != 0)
            {
                QString Name = Item->text();
                int n, N = S->getNumFitDataSets();
                for (n=0; n<N; ++n) if (S->getFitDataName(n) == Name)
                {
                    FitDataBox->addItem(Name);
                    break;
                }
            }
        }
        delete Item;
    }
}

void CoupledPotFitOutputImportDialog::OK()
{
     if (Mol != 0)
     {
        bool OK = true;
        if (Update->isChecked())
        {
            int n, m, N = ChoosenFitDataSets->count();
            for (n=m=0; n<N; ++n) if (Mol->getFitData(ChoosenFitDataSets->item(m)->text()) == 0)
            {
                OK = false;
                delete ChoosenFitDataSets->takeItem(m);
            }
        }
        if (OK) accept();
    }
    else
    {
        QMessageBox::information(this, "MolSpekAnalysis", "You have to select a molecule!");
    }
}
