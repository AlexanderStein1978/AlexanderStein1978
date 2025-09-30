//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "termview.h"
#include "molecule.h"
#include "atom.h"
#include "termtable.h"
#include "constants.h"
#include "utils.h"

#include <QStringList>

TermView::TermView(TermTable *Term) : TableWindow(TermEnergyView)
{
	int i, I1, I2, NSt = Term->getNumStates();
	QString N1, N2, Name = "Overview table for " + Term->getName();
	TT = Term;
	setName(Name);
	setSource(Term->getSource());
	Molecule *M = Term->getMolecule();
	Atom *A;
	if (M != 0)
	{
		A = M->getAtom1();
		if (A != 0) N1 = A->getChSymb();
		A = M->getAtom2();
		if (A != 0) N2 = A->getChSymb();
	}
	if (N1.isEmpty()) N1 = " ";
	Data = Term->getData();
	JMax = Term->getMaxJ();
	vMax = Term->getMaxv();
	NComp = Term->getNumComp();
	NIso = Term->getNumIso();
	View->addItems(QStringList() << "Term energies" << "Doublet separations" << "v differences"
				    << "Thermal populations");
	for (i=0; i < NSt; i++) View->addItem("Mix C" + QString::number(i));
	TUnit->addItem(QString("°C").right(2));
	TUnit->addItem("K");
	for (i=0; i < NIso; i++) 
	{
		Term->GetIsoZ(i, I1, I2);
		Iso->addItem(QString::number(I1) + N1 + QString::number(I2) + N2);
	}
	for (i=0; i < NComp; i++) Comp->addItem(QString::number(i));
	connect(View, SIGNAL(currentIndexChanged(int)), this, SLOT(ShowData()));
	connect(Iso, SIGNAL(currentIndexChanged(int)), this, SLOT(ShowData()));
	connect(Comp, SIGNAL(currentIndexChanged(int)), this, SLOT(ShowData()));
	connect(Temp, SIGNAL(textChanged(QString)), this, SLOT(ShowData()));
	connect(TUnit, SIGNAL(currentIndexChanged(int)), this, SLOT(ShowData()));
	ShowData();
}

TermView::~TermView()
{
}

void TermView::ShowData()
{
	int J, v, C = Comp->currentIndex(), I = Iso->currentIndex(), n, j;
	double D, T;
	Molecule *Mol = TT->getMolecule();
	ElState *St = TT->getElState();
	int JStart = (St != 0 ? St->getJStart(I, C) : 0), JStep = (Mol != 0 ? Mol->getJStep(I) : 1);
	QTableWidgetItem *It;
	switch(n = View->currentIndex())
	{
		case 0:
			Tab->setColumnCount(vMax + 1);
			Tab->setRowCount((JMax - JStart) / JStep + 1);
			for (v=0; v <= vMax; v++) 
				Tab->setHorizontalHeaderItem(v, new QTableWidgetItem("v=" + QString::number(v)));
			for (J = JStart, j=0; J <= JMax; J += JStep, j++)
			{
				Tab->setVerticalHeaderItem(j, new QTableWidgetItem("J=" + QString::number(J)));
				for (v=0; v<= vMax; v++) 
				{
					if ((It = Tab->item(j, v)) == 0) Tab->setItem(j, v, It = new QTableWidgetItem());
					It->setText(QString::number(Data[C][I][v][J], 'g', 11));
				}
			}
			Temp->setEnabled(false);
			TUnit->setEnabled(false);
			break;
		case 1:
			Tab->setColumnCount(vMax + 1);
			Tab->setRowCount((JMax - JStart) / JStep - 1);
			for (v=0; v <= vMax; v++)
				Tab->setHorizontalHeaderItem(v, new QTableWidgetItem("v''=" + QString::number(v)));
			for (J = JStart + 1, j=0; J < JMax; J += JStep, j++)
			{
				Tab->setVerticalHeaderItem(j, new QTableWidgetItem("J'=" + QString::number(J)));
				for (v=0; v <= vMax; v++)
				{
					if ((It = Tab->item(j, v)) == 0) 
						Tab->setItem(j, v, It = new QTableWidgetItem());
					It->setText(QString::number(
							(D = Data[C][I][v][J+1] - Data[C][I][v][J-1]) > 0.0 ? D : 0.0, 'g', 8));
				}
			}
			Temp->setEnabled(false);
			TUnit->setEnabled(false);
			break;
		case 2:
			Tab->setColumnCount(vMax);
			Tab->setRowCount((JMax - JStart) / JStep + 1);
			for (v=0; v < vMax; v++)
				Tab->setHorizontalHeaderItem(v, 
						new QTableWidgetItem("v=" + QString::number(v+1) + "-" + QString::number(v)));
			for (J = JStart, j=0; J <= JMax; J += JStep, j++)
			{
				Tab->setVerticalHeaderItem(j, new QTableWidgetItem("J=" + QString::number(J)));
				for (v=0; v < vMax; v++)
				{
					if ((It = Tab->item(j, v)) == 0) Tab->setItem(j, v, It = new QTableWidgetItem());
					It->setText(QString::number(
								(D = Data[C][I][v+1][J] - Data[C][I][v][J]) > 0.0 ? D : 0.0, 'g', 8));
				}
			}
			Temp->setEnabled(false);
			TUnit->setEnabled(false);
			break;
		case 3:
			Temp->setEnabled(true);
			TUnit->setEnabled(true);
			double **Pop;
			int nJ, nv;
			T = Temp->text().toDouble();
			if (TUnit->currentIndex() == 0) T += C_K_C;
			TT->getThermPopulation(Pop, nv, nJ, T, I, C);
			Tab->setColumnCount(vMax + 1);
			Tab->setRowCount((JMax - JStart) / JStep + 1);
			for (v=0; v <= vMax; v++) 
				Tab->setHorizontalHeaderItem(v, new QTableWidgetItem("v=" + QString::number(v)));
			for (J = JStart, j=0; J <= JMax; J += JStep, j++)
			{
				Tab->setVerticalHeaderItem(j, new QTableWidgetItem("J=" + QString::number(J)));
				for (v=0; v<= vMax; v++) 
				{
					if ((It = Tab->item(j, v)) == 0) Tab->setItem(j, v, It = new QTableWidgetItem());
					It->setText(QString::number(Pop[v][J], 'g', 5));
				}
			}
			Destroy(Pop, nv);
			break;
		default:
			Temp->setEnabled(false);
			TUnit->setEnabled(false);
			double *****MixCoeff = TT->getMixCoefficients();
			n-=4;
			Tab->setColumnCount(vMax + 1);
			Tab->setRowCount((JMax - JStart) / JStep + 1);
			for (v=0; v <= vMax; v++) Tab->setHorizontalHeaderItem(v, new QTableWidgetItem("v=" + QString::number(v)));
			for (J = JStart, j=0; J <= JMax; J += JStep, j++)
			{
				Tab->setVerticalHeaderItem(j, new QTableWidgetItem("J=" + QString::number(J)));
				for (v=0; v<= vMax; v++)
				{
					if ((It = Tab->item(j, v)) == 0) 
						Tab->setItem(j, v, new QTableWidgetItem(QString::number(MixCoeff[C][I][v][J][n], 'f', 5)));
					else It->setText(QString::number(MixCoeff[C][I][v][J][n], 'f', 5));
					
				}
			}
			break;
	}
}
