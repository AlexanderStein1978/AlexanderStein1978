//
// C++ Implementation: CouplingFuncs
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "CouplingFuncs.h"


CouplingFuncs::CouplingFuncs(MainWindow *MW) : Potential(MW)
{
    setType(MDIChild::CouplingFunction);
}

void CouplingFuncs::setData(int Nb, double *b, double binf, double Rs, double Rc, double epsilon, double xi_Rx_s, double xi_eps)
{
    m_Nb = Nb;
    Tab->blockSignals(true);
    Tab->setRowCount(Nb + 6);
    int n;
    for (n=0; n<Nb; ++n)
    {
        Tab->setItem(n, 0, new QTableWidgetItem("b_" + QString::number(n)));
        Tab->setItem(n, 1, new QTableWidgetItem(QString::number(b[n], 'g', 16)));
    }
    Tab->setItem(n, 0, new QTableWidgetItem("b_inf"));
    Tab->setItem(n, 1, new QTableWidgetItem(QString::number(binf, 'g', 16)));
    Tab->setItem(n+1, 0, new QTableWidgetItem("Rs"));
    Tab->setItem(n+1, 1, new QTableWidgetItem(QString::number(Rs)));
    Tab->setItem(n+2, 0, new QTableWidgetItem("Rc"));
    Tab->setItem(n+2, 1, new QTableWidgetItem(QString::number(Rc)));
    Tab->setItem(n+3, 0, new QTableWidgetItem("epsilon"));
    Tab->setItem(n+3, 1, new QTableWidgetItem(QString::number(epsilon)));
    Tab->setItem(n+4, 0, new QTableWidgetItem("x_Rx_s"));
    Tab->setItem(n+4, 1, new QTableWidgetItem(QString::number(xi_Rx_s)));
    Tab->setItem(n+5, 0, new QTableWidgetItem("xi_eps"));
    Tab->setItem(n+5, 1, new QTableWidgetItem(QString::number(xi_eps)));
    for (n=0; n <= Nb + 5; ++n)
    {
        Tab->setItem(n, 2, new QTableWidgetItem(""));
        Tab->setItem(n, 3, new QTableWidgetItem(""));
    }
    Tab->blockSignals(false);
    Changed();
}

void CouplingFuncs::getData(int &Nb, double *&b, double &binf, double &Rs, double &Rc, double &epsilon)
{
    Nb = m_Nb;
    b = new double[Nb];
    int n;
    for (n=0; n < Nb; ++n) b[n] = Tab->item(n, 1)->text().toDouble();
    binf = Tab->item(n, 1)->text().toDouble();
    Rs = Tab->item(n+1, 1)->text().toDouble();
    Rc = Tab->item(n+2, 1)->text().toDouble();
    epsilon = Tab->item(n+3, 1)->text().toDouble();
}

QTextEdit *CouplingFuncs::getTexTable()
{
    int Nb, NR = (m_Nb + 5) / 2, n;
    double *b, binf, Rs, Rc, epsilon;
    getData(Nb, b, binf, Rs, Rc, epsilon);
    QStringList TexTable(QStringList() << "\\begin{table}" << "\\caption{}\\\\" << "\\label{" + getName() + "}" << "\\centering" << "\\begin{tabular}{lr||lr}\\hline\\noalign{\\smallskip}");
    QStringList Data;
    QString Buffer;
    for (n=0; n < Nb; ++n)
    {
        Buffer = QString("$b_%1$ & %2").arg(n).arg(b[n], 0, 'f', 6);
        if (n>0) Buffer += QString(" \\AA$^{-%1}$").arg(n);
        Data << Buffer;
    }
    Data << QString("$b_\\infty$ & %1").arg(binf, 0, 'g', 6);
    Data << QString("$R_s$ & %1 \\AA").arg(Rs, 0, 'g', 6);
    Data << QString("$R_c$ & %1 \\AA").arg(Rc, 0, 'g', 6);
    Data << QString("$\\epsilon$ & %1").arg(epsilon, 0, 'g', 6);
    for (n = Nb; n < Data.size(); ++n) if (Data[n].indexOf('.') == -1) Data[n] += ".0";
    for (n=0; n < NR; ++n) TexTable << QString("%1 & %2 \\\\").arg(Data[n], (NR + n < Data.size() ? Data[NR + n] : "   "));
    TexTable << "\\hline" << "\\end{tabular}" << "\\end{table}";
    QTextEdit *Window = new QTextEdit();
    Window->setPlainText(TexTable.join("\n"));
    return Window;
}
