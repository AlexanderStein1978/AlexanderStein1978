//
// C++ Implementation: MeasuredTermEnergies
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "measuredtermenergies.h"


MeasuredTermEnergies::MeasuredTermEnergies(MainWindow *MW) : TableWindow(MDIChild::LineTab, MW)
{
	Data = 0;
	NumComp = NumIso = Maxv = MaxJ = numVRows = 0;
	Tab->setColumnCount(9);
	Tab->setHorizontalHeaderLabels(QStringList() << "Iso 1" << "Iso 2" << "component" 
			<< "v" << "J" << "term energy" << "Intensity" << "I rel. to sec." << "I rel. to av.");
}

void MeasuredTermEnergies::addData(int I, int Comp, int nTerm, double **nData, int **Ass)
{
	IsoTab *Iso = (molecule != 0 ? molecule->getIso() : 0);
	int n, in1, in2, r;
	bool iv = false, cv = false;
	if (Data != 0) 
	{
		Destroy (Data, NumComp, NumIso, Maxv + 1);
		Data = 0;
	}
	if (Iso != 0 ? I < Iso->numIso : false)
	{
		in1 = Iso->mNumIso1[I];
		in2 = Iso->mNumIso2[I];
	}
	else 
	{
		in1 = I;
		in2 = 0;
	}
	for (n=0; n < numVRows && (!iv || !cv); n++)
	{
		if (!iv) if (Tab->item(n, 0)->text().toInt() == in1 && Tab->item(n, 1)->text().toInt() == in2) 
			iv = true;
		if (!cv) if (Tab->item(n, 2)->text().toInt() == Comp) cv = true;
	}
	if (!iv) NumIso++;
	if (!cv) NumComp++;
	Tab->setRowCount(numVRows + nTerm);
	for (n=0, r = numVRows; n < nTerm; n++, r++)
	{
		Tab->setItem(r, 0, new QTableWidgetItem(QString::number(in1)));
		Tab->setItem(r, 1, new QTableWidgetItem(QString::number(in2)));
		Tab->setItem(r, 2, new QTableWidgetItem(QString::number(Comp)));
		Tab->setItem(r, 3, new QTableWidgetItem(QString::number(Ass[n][0])));
		Tab->setItem(r, 4, new QTableWidgetItem(QString::number(Ass[n][1])));
		Tab->setItem(r, 5, new QTableWidgetItem(QString::number(nData[n][0], 'g', 11)));
		Tab->setItem(r, 6, new QTableWidgetItem(QString::number(nData[n][1], 'g', 5)));
		Tab->setItem(r, 7, new QTableWidgetItem(QString::number(nData[n][2], 'g', 5)));
		Tab->setItem(r, 8, new QTableWidgetItem(QString::number(nData[n][3], 'g', 5)));
		if (Ass[n][0] > Maxv) Maxv = Ass[n][0];
		if (Ass[n][1] > MaxJ) MaxJ = Ass[n][1];
	}
	numVRows += nTerm;
	Changed();
}

double ****MeasuredTermEnergies::getData()
{
	if (Data == 0)
	{
		int C, I=0, v, J, n, I1, I2 = 0, Iso[NumIso][2], MI = 0;
		Data = Create(NumComp, NumIso, Maxv + 1, MaxJ + 1);
		for (C=0; C < NumComp; C++) for (I=0; I < NumIso; I++) for (v=0; v <= Maxv; v++)
					for (J=0; J <= MaxJ; J++) Data[C][I][v][J] = 0.0;
		for (n=0; n < numVRows; n++)
		{
			if ((I1 = Tab->item(n, 0)->text().toInt()) != Iso[I][0] 
						  || (I2 = Tab->item(n, 1)->text().toInt()) != Iso[I][1])
			{
				for (I=0; (I < MI ? I1 != Iso[I][0] || Iso[I][1] != I2 : false); I++) ;
				if (I < NumIso)
				{
					if (I == MI)
					{
						Iso[I][0] = I1;
						Iso[MI++][1] = I2;
					}
				}
				else
				{
					printf("MeasuredTermEnergies error: too many isotopologues!\n");
					continue;
				}
			}
			if ((v = Tab->item(n, 3)->text().toInt()) > Maxv)
			{
				printf("MeasuredTermEnergies error: too many vibrational levels!\n");
				continue;
			}
			if ((C = Tab->item(n, 2)->text().toInt()) >= NumComp)
			{
				printf("MeasuredTermEnergies error: too many components!\n");
				continue;
			}
			if ((J = Tab->item(n, 3)->text().toInt()) > MaxJ)
			{
				printf("MeasuredTermEnergies error: too many rotational levels!\n");
				continue;
			}
			Data[C][I][v][J] = Tab->item(n, 4)->text().toDouble();
		}
	}
	return Data;
}

bool MeasuredTermEnergies::readData(QString Filename)
{
	if (!TableWindow::readData(Filename)) return false;
	int Iso[MaxIso][2], I, I1, I2 = 0, r, R = Tab->rowCount(), b, C[10], c;
	NumIso = NumComp = 0;
	MaxJ = Maxv = -1;
	for (r=I=c=0; r<R; r++)
	{
		if ((I1 = Tab->item(r, 0)->text().toInt()) != Iso[I][0] 
				   || (I2 = Tab->item(r, 1)->text().toInt()) != Iso[I][1] || NumIso == 0)
		{
			for (I=0; (I < NumIso ? I1 != Iso[I][0] || I2 != Iso[I][1] : false); I++) ;
			if (I == MaxIso) 
			{
				printf("MeasuredTermEnergies error: too many isotopologues!\n");
				I--;
			}
			if (I == NumIso)
			{
				Iso[NumIso++][0] = I1;
				Iso[I][1] = I2;
			}
		}
		if ((b = Tab->item(r, 2)->text().toInt()) != C[c] || NumComp == 0)
		{
			for (c=0; (c < NumComp ? C[c] != b : false); c++) ;
			if (c == 10) 
			{
				printf("MeasuredTermEnergies error: too many components!\n");
				c--;
			}
			if (c == NumComp) C[NumComp++] = b;
		}
		if ((b = Tab->item(r, 3)->text().toInt()) > Maxv) Maxv = b;
		if ((b = Tab->item(r, 4)->text().toInt()) > MaxJ) MaxJ = b;
	}
	return true;
}

void MeasuredTermEnergies::setMolecule(Molecule *Mol)
{
	TableWindow::setMolecule(Mol);
}

bool MeasuredTermEnergies::writeData(QString FileName)
{
	return TableWindow::writeData(FileName);
}
