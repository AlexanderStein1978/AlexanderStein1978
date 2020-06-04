//
// C++ Implementation: dataplot
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2016
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "dataplot.h"
#include "MainWindow.h"
#include "linetable.h"
#include "elstate.h"
#include "termtable.h"
#include "utils.h"
#include "molecule.h"

#include <stdio.h>

#include <QComboBox>
#include <QPaintDevice>
#include <QPainter>
#include <QColor>
#include <QPen>
#include <QBrush>
#include <QPointF>
#include <QRect>


DataPlot::DataPlot(MainWindow *MW) : DiagWindow(DataSetPlot, MW)
{
	Block = false;
	Data = 0;
	NPoints = Nv = NJ = 0;
	state = 0;
	mol = 0;
	lineT = 0;
	termT = 0;
	XMin = -1.0;
	YMin = 0.0;
	setWindowTitle("Data field");
	xStartLabel->setText("v min:");
	xStopLabel->setText("v max:");
	yStartLabel->setText("J min:");
	yStopLabel->setText("J max:");
	XUnit = "v";
	YUnit = "J";
	SymbolSize = 5;
	SymbolCount = 3;
	SymbolColors = new QColor[SymbolCount];
	SymbolColors[0].setRgb(0, 0, 0);
	SymbolColors[1].setRgb(0, 0, 255);
	SymbolColors[2].setRgb(255, 0, 0);
	MolBox->setEditable(false);
	StateBox->setEditable(false);
	SourceBox->setEditable(false);
	connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(molBoxChanged()));
	connect(StateBox, SIGNAL(currentIndexChanged(int)), this, SLOT(stateBoxChanged(int)));
	connect(SourceBox, SIGNAL(currentIndexChanged(int)), this, SLOT(sourceBoxChanged()));
	moleculesChanged();
}


DataPlot::~DataPlot()
{
	if (Data != 0) Destroy(Data, NPoints);
	if (mol != 0) disconnect(mol, SLOT(propertiesChanged()), this, SLOT(molBoxChanged()));
	if (lineT != 0) disconnect(lineT, SLOT(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	if (termT != 0) disconnect(termT, SLOT(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	delete[] SymbolColors;
}

void DataPlot::drawSymbol(QPainter *P, int x, int y, int NR)
{
	int SD = (SymbolSize - 1) / 2;
	QPointF p[3];
	switch (NR)
	{
		case 1:
			P->setPen(SymbolColors[1]);
			P->setBrush(SymbolColors[1]);
			p[0] = QPointF(x, y - SD);
			p[1] = QPointF(x + SD, y + SD);
			p[2] = QPointF(x - SD, y + SD);
			P->drawConvexPolygon(p, 3);
			break;
		case 2:
			P->setPen(SymbolColors[2]);
			P->setBrush(SymbolColors[2]);
			p[0] = QPointF(x - SD, y - SD);
			p[1] = QPointF(x + SD, y - SD);
			p[2] = QPointF(x, y + SD);
			P->drawConvexPolygon(p, 3);
			break;
		default:
			P->fillRect(QRect(x - SD, y - SD, SymbolSize, SymbolSize), SymbolColors[0]);
			break;
	}
}

int DataPlot::heapSort(int sortFuncs(int *L1, int *L2), int **&LineData, int N)
{
	int i, n, m, L=0;
	int *B = LineData[0], **LB = new int*[N];
	while (B!=0)
	{
		B = 0;
		for (n=0, m=2; m<=N; n++, m+=2) 
		{
			if (sortFuncs(LineData[m-1], LineData[n]) == 1)
			{
				B = LineData[m-1];
				LineData[m-1] = LineData[n];
				LineData[n] = B;
			}
			if (m < N ? sortFuncs(LineData[m], LineData[n]) == 1 : false)
			{
				B = LineData[m];
				LineData[m] = LineData[n];
				LineData[n] = B;
			}
		}
	}
	for (i=0; i<N; i++)
	{
		if (L>0 ? sortFuncs(LB[L-1], LineData[0]) == 0 : false) delete[] LineData[0];
		else LB[L++] = LineData[0];
		for (n=0, m=2; LineData[n] != 0; m=2*((n=m)+1))
		{
			if (m < N)
			{
				if (sortFuncs(LineData[m-1], LineData[m]) == 1) m--;
			}
			else if (m==N) m--;
			else break;
			LineData[n] = LineData[m];
		}
		if (m >= N) LineData[n] = 0;
	}
	delete[] LineData;
	LineData = LB;
	return L;
}

void DataPlot::molBoxChanged()
{
	if (Block) return;
	Molecule *nmol;
	if (mol != (nmol = MW->getMolecule(MolBox->currentText())))
	{
		if (mol != 0) disconnect(mol, SLOT(propertiesChanged()), this, SLOT(molBoxChanged()));
		if ((mol = nmol) != 0) connect(mol, SLOT(propertiesChanged()), this, SLOT(molBoxChanged()));
	}
	Block = true;
	StateBox->clear();
	int m=0;
	if (mol != 0)
	{
		int N = mol->getNumStates(), n;
		ElState *S;
		for (n=0; n<N; n++)
		{
			StateBox->addItem((S = mol->getStateP(n))->getName());
			if (S == state) StateBox->setCurrentIndex(m=n);
		}
	}
	Block = false;
	stateBoxChanged(m);
}

void DataPlot::moleculesChanged()
{
	Block = true;
	MolBox->clear();
	int n, m=-1, N = MW->getNumMolecules();
	Molecule *mb;
	for (n=0; n<N; n++) 
	{
		MolBox->addItem((mb = MW->getMolecule(n))->getName());
		if (mb == mol) MolBox->setCurrentIndex(m=n);
	}
	Block = false;
	if (m == -1) molBoxChanged();
}

void DataPlot::PSpektrum(QPainter &P, const QRect &A, bool /*PrintFN*/)
{
	if (Block) return;
	int SD = (SymbolSize - 1) / 2, n, x, y;
	for (n=0; n < NPoints; n++) 
		if ((x = int(XO + XSF * double(Data[n][1]))) > ScaleYWidth + SD
				   && (y = int(YO - YSF * double(Data[n][2]))) < A.height() - ScaleXHeight - SD)
			drawSymbol(&P, x, y, Data[n][0]);
}

void DataPlot::sourceBoxChanged()
{
	if (Block) return;
	int n, l, N, T, I, v, J, C, NI, NC, Mv=-1, MJ=-1;
	double ****TData;
	QString Name = SourceBox->currentText();
	TermTable *TT = 0;
	LineTable *LT = 0;
	Transition *Tr = 0;
	if (Data != 0) Destroy(Data, NPoints);
	Nv = NJ = 0;
	NPoints = 0;
	if (state == 0) Data = 0;
	else 
	{
		for (n=0, N = state->getNumTermTables(); (n<N ? 
				   (TT = state->getTermTable(n))->getName() != Name : false); n++) ;
		if (n<N)
		{
			TData = TT->getData();
			NC = TT->getNumComp();
			NI = TT->getNumIso();
			Nv = TT->getMaxv() + 1;
			NJ = TT->getMaxJ() + 1;
			for (I=0; I < NI; I++) for (v=0; v < Nv; v++) for (J=0; J < NJ; J++)
			{
				for (C=0; (C < NC ? TData[C][I][v][J] == 0.0 : false); C++) ;
				if (C < NC) NPoints++;
			}
			Data = CreateInt(NPoints, 3);
			for (I = NPoints = 0; I < NI; I++) for (v=0; v < Nv; v++) for (J=0; J < NJ; J++)
			{
				for (C=0; (C < NC ? TData[C][I][v][J] == 0.0 : false); C++) ;
				if (C < NC)
				{
					if (v > Mv) Mv = v;
					if (J > MJ) MJ = J;
					Data[NPoints][0] = I;
					Data[NPoints][1] = v;
					Data[NPoints++][2] = J;
				}
			}
		}
		else
		{
			TT = 0;
			for (n=l=N=0, T = mol->getNumTransitions(); n<T && l==N; n++)
				for (l=0, N = (Tr = mol->getTransitionP(n))->getNumLineTables(); 
								 (l<N ? Tr->getLineTableName(l) != Name : false); l++) ;
			if (l<N)
			{
				LT = Tr->getLineTable(l);
				double *E = new double[N = LT->getAnzahlLinien()];
				int **LineData = CreateInt(N, 6);
				double *Unc = new double[N];
				LT->getLines(LineData, E, Unc);
				delete[] E;
				delete[] Unc;
				//printf("Vor heapsort\n");
				if (Tr->getLowerState() == state)
				{
					N = heapSort(sortLS, LineData, N);
					l = 3;
				}
				else
				{
					N = heapSort(sortUS, LineData, N);
					l = 1;
				}
				//printf("Nach heapsort\n");
				Data = CreateInt(N, 3);
				for (n=0; n<N; n++) 
				{
					Data[n][0] = LineData[n][0];
					Data[n][1] = LineData[n][l];
					Data[n][2] = LineData[n][l+1];
					if (Data[n][1] > Mv) Mv = Data[n][1];
					if (Data[n][2] > MJ) MJ = Data[n][2];
				}
				Destroy(LineData, N);
				//printf("Nach destroy\n");
			}
			else 
			{
				LT = 0;
				Data = 0;
			}
		}
	}
	if (lineT != 0) disconnect(lineT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	else 
		if (termT != 0) disconnect(termT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	if (LT != 0) connect(lineT = LT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	else if (TT != 0) connect(termT = TT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
	XMax = 1.5 * double(Mv);
	YMax = 1.5 * double(MJ);
	xStart->setText("-1");
	xStop->setText(QString::number(1.1 * double(Mv)));
	yStart->setText("0");
	yStop->setText(QString::number(1.1 * double(MJ)));
	printf("Vor paint\n");
	Paint();
}

void DataPlot::stateBoxChanged(int i)
{
	if (Block) return;
	if (mol == 0)
	{
		state = 0;
		SourceBox->clear();
		return;
	}
	ElState *St = mol->getStateP(i);
	Transition *Tr;
	if (St == state) return;
	Block = true;
	int T = mol->getNumTransitions(), n, m, N;
	SourceBox->clear();
	for (n=0, N = St->getNumTermTables(); n<N; n++) SourceBox->addItem(St->getTermTableName(n));
	for (n=0; n<T; n++) 
	{
		Tr = mol->getTransitionP(n);
		if (Tr->getUpperState() == St || Tr->getLowerState() == St)
			for (m=0, N = Tr->getNumLineTables(); m<N; m++) 
				SourceBox->addItem(Tr->getLineTableName(m));
	}
	state = St;
	Block = false;
	sourceBoxChanged();
}
