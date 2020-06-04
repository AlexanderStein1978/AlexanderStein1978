//
// C++ Implementation: PerturbationInfoTable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "PerturbationInfoTable.h"


PerturbationInfoTable::PerturbationInfoTable(QWidget *i_parent) : QWidget(i_parent), m_Tab(new PerturbationInfoTab(this)),
    m_PerturbationSelectBox(new QComboBox(this)), m_ResidualFit(0)
{
    setWindowTitle("Properties of local perturbation");
    QGridLayout* Layout = new QGridLayout(this);
    Layout->addWidget(new QLabel("Selected:", this), 0, 0);
    Layout->addWidget(m_PerturbationSelectBox, 0, 1);
    Layout->setRowStretch(0, 0);
    Layout->addWidget(m_Tab, 1, 0, 1, 2);
    Layout->setRowStretch(1, 100);
    m_PerturbationSelectBox->setEditable(false);
    connect(m_PerturbationSelectBox, SIGNAL(currentIndexChanged(int)), this, SLOT(SelectPerturbation(int)));
    m_Tab->setColumnCount(3);
    m_Tab->setHorizontalHeaderLabels(QStringList() << "Parameter:" << "Value:" << "Uncertainty:");
    connect(m_Tab, SIGNAL(RightClicked(QPoint)), this, SLOT(ShowPopupMenu(QPoint)));
}

PerturbationInfoTable::~PerturbationInfoTable()
{
}

void PerturbationInfoTable::SelectPerturbation(int i_index)
{
    if (i_index < 0) return;
    if (i_index != m_PerturbationSelectBox->currentIndex())
    {
        m_PerturbationSelectBox->setCurrentIndex(i_index);
        return;
    }
    const LocalPerturbation* currentPerturbation = m_ResidualFit->getLocalPerturbation(i_index);
    int NBandConstants = currentPerturbation->GetNBandC();
    double Value, Uncertainty;
    m_Tab->setRowCount(NBandConstants + 3);
    m_Tab->setItem(0, 0, new QTableWidgetItem("Center J:"));
    m_Tab->setItem(0, 1, new QTableWidgetItem(QString::number(currentPerturbation->GetCenter(), 'f', 1)));
    m_Tab->setItem(0, 2, new QTableWidgetItem("0.5"));
    m_Tab->setItem(1, 0, new QTableWidgetItem("H_12"));
    currentPerturbation->GetH12(Value, Uncertainty);
    m_Tab->setItem(1, 1, new QTableWidgetItem(QString::number(Value, 'g', 6)));
    m_Tab->setItem(1, 2, new QTableWidgetItem(QString::number(Uncertainty, 'g', 6)));
    for (int i=0; i < NBandConstants; ++i)
    {
        currentPerturbation->GetBandC(i, Value, Uncertainty);
        m_Tab->setItem(i+2, 0, new QTableWidgetItem("Y_x" + QString::number(i)));
        m_Tab->setItem(i+2, 1, new QTableWidgetItem(QString::number(Value, 'g', 6)));
        m_Tab->setItem(i+2, 2, new QTableWidgetItem(QString::number(Uncertainty, 'g', 6)));
    }
    m_Tab->setItem(NBandConstants + 2, 0, new QTableWidgetItem("Sigma"));
    m_Tab->setItem(NBandConstants + 2, 1, new QTableWidgetItem(QString::number(currentPerturbation->GetFitSigma(), 'g', 6)));
    m_Tab->setItem(NBandConstants + 2, 2, new QTableWidgetItem);
    for (int r=0; r < m_Tab->rowCount(); ++r) for (int c=0; c < 3; ++c) m_Tab->item(r, c)->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
}

void PerturbationInfoTable::SetResidualFit(ResidualFit * const i_ResidualFit)
{
    m_ResidualFit = i_ResidualFit;
    int n, NPerturbations = (m_ResidualFit != 0 ? m_ResidualFit->getNumberOfLocalPerturbations() : 0);
    while ((n = m_PerturbationSelectBox->count()) < NPerturbations) m_PerturbationSelectBox->addItem("perturbation " + QString::number(n));
    while ((n = m_PerturbationSelectBox->count()) > NPerturbations) m_PerturbationSelectBox->removeItem(n-1);
    m_PerturbationSelectBox->setEnabled(NPerturbations > 1);
}

void PerturbationInfoTable::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton && m_ResidualFit != 0)
    {
        int currentPerturbationIndex = (m_ResidualFit->getNumberOfLocalPerturbations() > 1 ? m_PerturbationSelectBox->currentIndex() : 0);
        const LocalPerturbation* perturbation = m_ResidualFit->getLocalPerturbation(currentPerturbationIndex);
        if (0 != perturbation && perturbation->GetNBandC() > 2)
        {
            event->accept();
            ShowPopupMenu(event->pos());
        }
    }
    else
    {
        QWidget::mouseReleaseEvent(event);
    }
}

void PerturbationInfoTable::ShowPopupMenu(QPoint i_pos)
{
    QMenu popupMenu;
    QAction *reduceCoefficientsAction = new QAction("&Reduce band constants...", this);
    reduceCoefficientsAction->setStatusTip("Reduce the count of band constants to a number to select.");
    popupMenu.addAction(reduceCoefficientsAction);
    connect(reduceCoefficientsAction, SIGNAL(triggered()), this, SLOT(ReduceNumberOfBandConstants()));
    popupMenu.exec(i_pos);
}

void PerturbationInfoTable::ReduceNumberOfBandConstants()
{
   int currentPerturbationIndex = (m_ResidualFit->getNumberOfLocalPerturbations() > 1 ? m_PerturbationSelectBox->currentIndex() : 0);
   LocalPerturbation* perturbation = m_ResidualFit->getLocalPerturbation(currentPerturbationIndex);
   int N = perturbation->GetNBandC();
   QDialog* dialog = new QDialog(this);
   QGridLayout* layout = new QGridLayout(dialog);
   layout->addWidget(new QLabel("New number of band constants:", dialog), 0, 0);
   QComboBox *box = new QComboBox(dialog);
   for (int n=2; n<N; ++n) box->addItem(QString::number(n));
   box->setEditable(false);
   if (N==3) box->setEnabled(false);
   layout->addWidget(box, 0, 1);
   layout->setRowMinimumHeight(1, 20);
   QPushButton* okButton = new QPushButton("OK", dialog), *cancelButton = new QPushButton("Cancel", dialog);
   layout->addWidget(okButton, 2, 0);
   layout->addWidget(cancelButton, 2, 1);
   connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
   connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
   if (dialog->exec() == QDialog::Accepted)
   {
       double pp = perturbation->GetCenter(), Ec = perturbation->GetPointUpOBCa(pp);
       int n = box->currentIndex() + 2;
       for (int i=N; i>n; i--) perturbation->RemoveLastBandC();
       double En = perturbation->GetPointCalcBand(pp), T, Buff;
       perturbation->GetBandC(0, T, Buff);
       perturbation->SetBandC(0, T + Ec - En, Buff);
       perturbation->LevenbergMarquardt(100, 1e-6);
       m_ResidualFit->DetermineUncertaintiesByBeVariation(*perturbation);
       SelectPerturbation(currentPerturbationIndex);
   }
   delete dialog;
}
