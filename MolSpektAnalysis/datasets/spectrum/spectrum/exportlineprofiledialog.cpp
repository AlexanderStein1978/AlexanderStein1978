//
// C++ Implementation: spektrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "exportlineprofiledialog.h"


ExportLineProfileDialog::ExportLineProfileDialog(Spektrum * const i_spectrum, MainWindow* i_MW) : QWidget(i_MW) , m_spectrum(i_spectrum)
{
    QGridLayout* L = new QGridLayout(this);
    L->addWidget(new QLabel("Select line E_c [cm^-1]:", this), 0, 0);
    L->addWidget(m_lineSelectBox = new QComboBox(this), 0, 1, 1, 2);
    L->addWidget(new QLabel("File name:", this), 1, 0);
    L->addWidget(m_FileNameEdit = new QLineEdit("LineProfile.dat", this), 1, 1);
    L->addWidget(m_selectFileButton = new QPushButton("...", this), 1, 2);
    L->addWidget(new QLabel("Min. rel. E [cm^-1]:"), 2, 0);
    L->addWidget(m_eStartEdit = new QLineEdit("0.0", this), 2, 1, 1, 2);
    L->addWidget(new QLabel("Max. rel. E [cm^-1]:"), 3, 0);
    L->addWidget(m_eEndEdit = new QLineEdit("2.0", this), 3, 1, 1, 2);
    L->addWidget(new QLabel("Step size [cm^-1]:", this), 4, 0);
    L->addWidget(m_eStepEdit = new QLineEdit("1e-5", this), 4, 1, 1, 2);
    L->setRowMinimumHeight(5, 20);
    L->addWidget(m_OKButton = new QPushButton("OK", this), 6, 0);
    L->addWidget(m_CancelButton = new QPushButton("Cancel", this), 6, 1, 1, 2);

    int N = m_spectrum->GetNumFittedLines();
    for (int n = 0; n < N; ++n) m_lineSelectBox->addItem(QString::number(m_spectrum->GetFittedLine(n)->GetEcenter(), 'f', 4));
    m_lineSelectBox->setEditable(false);
    m_eStartEdit->setValidator(new QDoubleValidator);
    m_eEndEdit->setValidator(new QDoubleValidator);
    m_eStepEdit->setValidator(new QDoubleValidator);

    connect(m_selectFileButton, SIGNAL(clicked()), this, SLOT(SelectFileName()));
    connect(m_lineSelectBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LineChanged(int)));
    connect(m_OKButton, SIGNAL(clicked()), this, SLOT(Write()));
    connect(m_CancelButton, SIGNAL(clicked()), this, SIGNAL(closeThis()));
}

void ExportLineProfileDialog::SelectFileName()
{
    QString newFileName = QFileDialog::getOpenFileName(this, "Please select file name", m_FileNameEdit->text());
    if (!newFileName.isEmpty()) m_FileNameEdit->setText(newFileName);
}

void ExportLineProfileDialog::LineChanged(int i_index)
{
    if (i_index < 0) return;
    const Gaussian* currentLine = m_spectrum->GetFittedLine(i_index);
    double Estart, Eend, B, E, Width, Offset, Istart, Iend;
    currentLine->GetValues(B, E, Width, Offset);
    currentLine->GetDataRange(Estart, Eend);
    Width = Eend - Estart;
    if (B > 0)
    {
        Istart = Offset - 0.1 * B;
        Iend = Offset + 1.1 * B;
    }
    else
    {
        Istart = Offset + 1.1 * B;
        Iend = Offset - 0.1 * B;
    }
    m_spectrum->setRanges(Estart - 0.1 * Width, Eend + 0.1 * Width, Istart, Iend);
}

void ExportLineProfileDialog::Write()
{
    int index = m_lineSelectBox->currentIndex();
    if (index < 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Invalid input values: No line selected.");
        return;
    }
    double Estart = m_eStartEdit->text().toDouble(), Eend = m_eEndEdit->text().toDouble(), Estep = m_eStepEdit->text().toDouble(), E;
    if (abs(Estep) < 1e-9 * abs(Eend - Estart))
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Invalid input values: Step size to small compared to E range.");
        return;
    }
    QFile File(m_FileNameEdit->text());
    if (!File.open(QIODevice::WriteOnly))
    {
        QMessageBox::information(this, "MolSpektAnalysis", "Could not open file for writing.");
        return;
    }
    QTextStream S(&File);
    bool backward = Estart > Eend || Estep < 0.0;
    if (Estart > Eend)
    {
        E = Estart;
        Estart = Eend;
        Eend = E;
    }
    if (Estep < 0.0) Estep *= -1.0;
    int nEDig = static_cast<int>(-log10(Estep)) + 1;
    if (nEDig < 4) nEDig = 4;
    const Gaussian* L = m_spectrum->GetFittedLine(index);
    S << "E [cm^-1]\tI [a.u.]\n";
    if (backward) for (E = Eend; E >= Estart; E -= Estep) S << QString("%1\t%2\n").arg(E, 0, 'f', nEDig).arg(L->GetProfilePoint(E), 0, 'f', 4);
    else for (E = Estart; E <= Eend; E+= Estep) S << QString("%1\t%2\n").arg(E, 0, 'f', nEDig).arg(L->GetProfilePoint(E), 0, 'f', 4);
    emit closeThis();
}
