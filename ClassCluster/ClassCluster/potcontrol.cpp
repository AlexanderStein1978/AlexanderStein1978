#include "potcontrol.h"
#include "potstruct.h"
#include "potential.h"
#include "potentialplot.h"
#include "controlwindow.h"

#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QDoubleValidator>
#include <QLabel>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>


PotControl::PotControl(ControlWindow *i_parent)
    : parent(i_parent)
    , pot(nullptr)
    , plot(nullptr)
    , fileName(new QLineEdit(parent))
    , VScale(new QLineEdit("1.0", parent))
    , RScale(new QLineEdit("1.0", parent))
    , openB(new QPushButton("...", parent))
    , saveB(new QPushButton("save", parent))
    , saveAsB(new QPushButton("save as...", parent))
    , adjustReB(new QPushButton("adjust Re", parent))
    , showBox(new QCheckBox("plot", parent))
    , changed(false)
    , changing(false)
{
    VScale->setValidator(new QDoubleValidator(1e-5, 1e5, 1000, VScale));
    RScale->setValidator(new QDoubleValidator(1e-5, 1e5, 1000, RScale));

    connect(fileName, SIGNAL(editingFinished()), this, SLOT(Open()));
    connect(VScale, SIGNAL(editingFinished()), this, SLOT(Changed()));
    connect(RScale, SIGNAL(editingFinished()), this, SLOT(Changed()));
    connect(openB, SIGNAL(clicked()), this, SLOT(ShowOpenDialog()));
    connect(saveB, SIGNAL(clicked()), this, SLOT(Save()));
    connect(saveAsB, SIGNAL(clicked()), this, SLOT(SaveAs()));
    connect(adjustReB, SIGNAL(clicked()), this, SLOT(adjustRe()));
    connect(showBox, SIGNAL(toggled(bool)), this, SLOT(Plot(bool)));
}

PotControl::~PotControl()
{
    if (nullptr != pot) delete pot;
}

void PotControl::Init(const QString& data)
{
    if (data.isEmpty()) return;
    QStringList list = data.split('\t');
    fileName->setText(list[0]);
    if (list.size() > 1) VScale->setText(list[1]);
    if (list.size() > 2) RScale->setText(list[2]);
    Open();
}

void PotControl::Serialize(QTextStream& stream) const
{
    stream << fileName->text() << '\t' << VScale->text() << '\t' << RScale->text() << '\n';
}

void PotControl::FillLayout(QGridLayout* layout, const int row) const
{
    layout->addWidget(fileName, row, 1);
    layout->addWidget(openB, row, 2);
    layout->addWidget(saveB, row, 3);
    layout->addWidget(saveAsB, row, 4);
    layout->addWidget(new QLabel("VS:", parent), row, 5);
    layout->addWidget(VScale, row, 6);
    layout->addWidget(new QLabel("RS:", parent), row, 7);
    layout->addWidget(RScale, row, 8);
    layout->addWidget(adjustReB, row, 9);
    layout->addWidget(showBox, row, 10);
}

void PotControl::FillStruct(PotStruct& potStruct) const
{
    potStruct.pot = pot;
    potStruct.VZoom = VScale->text().toDouble();
    potStruct.RZoom = RScale->text().toDouble();
}

bool PotControl::canPotBeClosed() const
{
    if (nullptr == pot || pot->isSaved()) return true;

    QMessageBox::StandardButton result = QMessageBox::question(parent, "The potential " + fileName->text() + " is not saved!", "Do you want to save it now?",
                                                               QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel, QMessageBox::Yes);
    switch(result)
    {
    case QMessageBox::Yes:
        if (!pot->writeData()) return false;
    case QMessageBox::No:
        return true;
    case QMessageBox::Cancel:
        return false;
    default:
        break;
    }
    return true;
}

void PotControl::Open()
{
    if (fileNBackUp != fileName->text() && canPotBeClosed())
    {
        if (fileName->text().isEmpty())
        {
            pot->close();
            pot->deleteLater();
            pot = nullptr;
        }
        else openPotential();
    }
}

void PotControl::ShowOpenDialog()
{
    if (!canPotBeClosed()) return;
    QString startPath = (fileName->text().isEmpty() ? "../../../Physics/ClassCluster/" : fileName->text());
    QString newPath = QFileDialog::getOpenFileName(parent, "Please select the potential file", startPath, "*.pot");
    if (newPath.isEmpty()) return;
    fileName->setText(newPath);
    openPotential();
}

void PotControl::openPotential()
{
    Potential* newPot = new Potential;
    if (newPot->readData(fileName->text()))
    {
        if (nullptr != plot && showBox->isChecked())
        {
            Plot(false);
            plot->addPotential(newPot);
        }
        if (pot != nullptr)
        {
            disconnect(pot, SIGNAL(propertiesChanged()), this, SLOT(RecalcExtensions()));
            pot->close();
            pot->deleteLater();
        }
        pot = newPot;
        connect(pot, SIGNAL(propertiesChanged()), this, SLOT(RecalcExtensions()));
        fileNBackUp = fileName->text();
        int n = fileNBackUp.lastIndexOf(DIRSEP) + 1;
        pot->setName(fileNBackUp.mid(n, fileNBackUp.indexOf('.', n) - n));
        pot->Saved();
        pot->show();
    }
    else
    {
        delete newPot;
        if (showBox->isChecked()) showBox->setChecked(false);
    }
    if (nullptr != pot) setRelativePath();
}

void PotControl::setRelativePath()
{
    QString path = pot->getRelativePath(parent->getProgramPath());
    if (path != fileNBackUp)
    {
        if (path.left(2) == "..") path = "../" + path;
        setFileName(path);
    }
}

void PotControl::Save()
{
    if (nullptr != pot) pot->writeData(fileName->text());
}

void PotControl::SaveAs()
{
    if (nullptr != pot)
    {
        QString startPath = (fileName->text().isEmpty() ? "../../../Physics/ClassCluster/" : fileName->text());
        QString newFileName = QFileDialog::getSaveFileName(parent, "Please select the potential file", startPath, "*.pot");
        if (!newFileName.isEmpty())
        {
            pot->writeData(newFileName);
            setRelativePath();
        }
    }
}

void PotControl::Plot(const bool show)
{
    if (nullptr == plot) return;
    if (show) plot->plotPotential(pot);
    else plot->removePotential(pot);
    if (!plot->isVisible())
    {
        plot->setShowHistory(false);
        plot->show();
    }
}

void PotControl::RecalcExtensions()
{
    if (nullptr == pot || changing) return;
    changing = true;
    pot->cdConnectSR();
    pot->cdConnectLR1C();
    changing = false;
}

void PotControl::adjustRe()
{
    double Re, De, tRe = parent->getRe();
    pot->getReDe(Re, De);
    RScale->setText(QString::number(Re / tRe, 'f'));
}

void PotControl::closePot()
{
    pot->close();
}

void PotControl::setFileName(const QString& newFileName)
{
    fileNBackUp = newFileName;
    fileName->setText(newFileName);
}
