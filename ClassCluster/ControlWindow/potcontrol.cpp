#include "potcontrol.h"
#include "potstruct.h"
#include "potential.h"
#include "potentialplot.h"
#include "controlwindow.h"
#include "MainWindow.h"

#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QDoubleValidator>
#include <QLabel>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>


PotControl::PotControl(ControlWindow *i_parent, MainWindow *mw)
    : parent(i_parent)
    , pot(nullptr)
    , VScale(new QLineEdit("1.0", parent))
    , RScale(new QLineEdit("1.0", parent))
    , PotentialBox(new QComboBox)
    , adjustReB(new QPushButton("adjust Re", parent))
    , showBox(new QCheckBox("plot", parent))
    , MW(mw)
    , changed(true)
    , changing(false)
{
    VScale->setValidator(new QDoubleValidator(1e-5, 1e5, 1000, VScale));
    RScale->setValidator(new QDoubleValidator(1e-5, 1e5, 1000, RScale));
    PotentialBox->setEditable(false);

    connect(VScale, SIGNAL(editingFinished()), this, SLOT(Changed()));
    connect(RScale, SIGNAL(editingFinished()), this, SLOT(Changed()));
    connect(adjustReB, SIGNAL(clicked()), this, SLOT(adjustRe()));
    connect(showBox, SIGNAL(toggled(bool)), this, SLOT(Plot(bool)));
    connect(PotentialBox, SIGNAL(currentIndexChanged(int)), this, SLOT(PotentialBoxIndexChanged(int)));

    UpdatePotentialBox();
}

PotControl::~PotControl()
{
}

void PotControl::Init(const QString& data)
{
    if (data.isEmpty()) return;
    QStringList list = data.split('\t');
    exchangePotential(MW->getPotential(list[0], nullptr));
    if (list.size() > 1) VScale->setText(list[1]);
    if (list.size() > 2) RScale->setText(list[2]);
}

bool PotControl::isChangedSinceLastRun()
{
    if (changed)
    {
        changed = false;
        return true;
    }
    return false;
}

void PotControl::Serialize(QTextStream& stream, const QString &programPath)
{
    stream << (pot != nullptr ? pot->getRelativePath(programPath + "/DummyString.dat") : "") << '\t' << VScale->text() << '\t' << RScale->text() << '\n';
}

void PotControl::FillLayout(QGridLayout* layout, const int row) const
{
    layout->addWidget(new QLabel("Potential:", parent), row, 1);
    layout->addWidget(PotentialBox, row, 2);
    layout->addWidget(new QLabel("VS:", parent), row, 3);
    layout->addWidget(VScale, row, 4);
    layout->addWidget(new QLabel("RS:", parent), row, 5);
    layout->addWidget(RScale, row, 6);
    layout->addWidget(adjustReB, row, 7);
    layout->addWidget(showBox, row, 8);
}

void PotControl::FillStruct(PotStruct& potStruct) const
{
    potStruct.pot = pot;
    potStruct.VZoom = VScale->text().toDouble();
    potStruct.RZoom = RScale->text().toDouble();
}

void PotControl::exchangePotential(Potential *const newPot)
{
    if (pot != nullptr) disconnect(pot, SIGNAL(propertiesChanged()), this, SLOT(RecalcExtensions()));
    if (showBox->isChecked())
    {
        Plot(false);
        pot = newPot;
        Plot(true);
    }
    else pot = newPot;
    connect(pot, SIGNAL(propertiesChanged()), this, SLOT(RecalcExtensions()));
}

void PotControl::PotentialBoxIndexChanged(const int newIndex)
{
    exchangePotential(MW->getPotential(newIndex));
}

void PotControl::Plot(const bool show)
{
    parent->showPotential(pot, show);
}

void PotControl::PLotCloses()
{
    showBox->blockSignals(true);
    showBox->setChecked(false);
    showBox->blockSignals(false);
}

void PotControl::RecalcExtensions()
{
    if (nullptr == pot || changing) return;
    changing = true;
    pot->cdConnectSR();
    pot->cdConnectLR1C();
    changing = false;
    changed = true;
}

void PotControl::UpdatePotentialBox()
{
    PotentialBox->blockSignals(true);
    PotentialBox->clear();
    for (int i=0; i < MW->getNumPotentials(); ++i)
    {
        Potential* curPot = MW->getPotential(i);
        PotentialBox->addItem(curPot->getName());
        if (pot == curPot) PotentialBox->setCurrentIndex(i);
    }
    PotentialBox->blockSignals(false);
    if (nullptr == pot) exchangePotential(MW->getPotential(0));
}

void PotControl::adjustRe()
{
    double Re, De, tRe = parent->getRe();
    pot->getReDe(Re, De);
    RScale->setText(QString::number(tRe / Re, 'f'));
}
