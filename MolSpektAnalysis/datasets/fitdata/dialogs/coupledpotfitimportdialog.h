//
// C++ Interface: CoupledPotFitOutputImportDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef COUPLEDPOTFITOUTPUTIMPORTDIALOG_H
#define COUPLEDPOTFITOUTPUTIMPORTDIALOG_H


#include <QDialog>


class QPushButton;
class QListWidget;
class QComboBox;
class QCheckBox;
class QLineEdit;

class Molecule;
class MainWindow;
class FitData;


class CoupledPotFitOutputImportDialog : public QDialog
{
    Q_OBJECT

public:
    CoupledPotFitOutputImportDialog(MainWindow *MW);
    void getResult(Molecule *&Mol, FitData **&fDats, int &NfDats, double &Tolerancy);

private slots:
    void MoleculeChanged(int Index);
    void UpdateChanged(bool update);
    void StateChanged(int Index);
    void Add();
    void Remove();
    void OK();

private:
    QPushButton *OKButton, *CancelButton, *AddButton, *RemoveButton;
    QListWidget *ChoosenFitDataSets;
    QComboBox *MolBox, *StateBox, *FitDataBox;
    QCheckBox *Update;
    QLineEdit *TolerancyEdit;

    MainWindow *MW;
    Molecule *Mol;
};

#endif
