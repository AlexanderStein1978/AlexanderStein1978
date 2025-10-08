//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "wavefuncplot.h"
#include "MainWindow.h"
#include "constants.h"
#include "potential.h"
#include "molecule.h"
#include "utils.h"
#include "elstate.h"
#include "CoupledSineWaveFunc.h"
#include "isotab.h"
#include "datensatz.h"

#include <stdio.h>
#include <math.h>

#include <QColor>
#include <QComboBox>
#include <QPaintDevice>
#include <QPainter>
#include <QLineEdit>

WaveFuncPlot::WaveFuncPlot(int numWFPoints, MainWindow *MW, Potential *pot, int Iso, int J, int v) : DiagWindow(MDIChild::WaveFunctionPlot, MW)
{
	//printf("WaveFuncPlot::WaveFuncPlot\n");
    NumWFPoints = numWFPoints;
    nDatenS = 1;
	NChannels = 0;
	Pot = new Potential*[10];
	E = 0;
	F = 0;
	N = 0;
	Nv = 0;
	ri = 0;
	ra = 0;
	P = 0;
	WF = 0;
	Ra = 0.0;
	Ro = 0.0;
	Mol = 0;
	PMin = PMax = 0.0;
	h = (rmax - rmin) / (NumWFPoints - 1);
	WaveFuncSF = 0;
	PotC = new QColor(0, 255, 0);
	RotPotC = new QColor(255, 0, 0);
	LevelC = new QColor(150, 150, 150);
	WaveC = new QColor(0, 0, 0);
	SWF = 0;
	XUnit = "R [A]";
	YUnit = "Energy [cm^{-1}]";
	//printf("Vor SpektColor\n");
	SpektColor = *WaveC;
	//printf("Nach SpektColor\n");
	int n, NP, NM = MW->getNumMolecules(), NS, p, s;
	Molecule *M;
	ElState *S;
	QString PN;
	for (n=0; n < NM; n++) for (s=0, NS = (M = MW->getMolecule(n))->getNumStates(); s < NS; s++)
			for (p=0, NP = (S = M->getStateP(s))->getNumPotentials(); p < NP; p++) 
				if (!(PN = S->getPotentialName(p)).isEmpty())
	{
		StateBox->addItem(PN);
		if (pot != 0 ? PN == pot->getName() : false) StateBox->setCurrentIndex(n);
	}
	if (Iso >= 0)
	{
		for (n=0; n <= Iso; n++) IsoBox->addItem(QString::number(n));
		IsoBox->setCurrentIndex(Iso);
		if (pot != 0) Mol = pot->getMolecule();
	}
	JBox->addItem(QString::number(J));
	if (v>=0)
	{
		for (n=0; n<=v; n++) vBox->addItem(QString::number(n));
		vBox->setCurrentIndex(v);
	}
	PotChanged(StateBox->currentText());
	connect(StateBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(PotChanged(QString)));
	connect(IsoBox, SIGNAL(currentIndexChanged(int)), this, SLOT(IsoChanged(int)));
	connect(JBox, SIGNAL(editTextChanged(QString)), this, SLOT(JChanged()));
	connect(JBox, SIGNAL(currentIndexChanged(int)), this, SLOT(JChanged()));
	connect(vBox, SIGNAL(currentIndexChanged(int)), this, SLOT(vChanged(int)));
	//printf("Ende WaveFuncPlot\n");
}


WaveFuncPlot::~WaveFuncPlot()
{
	//printf("WaveFuncPlot::~WaveFuncPlot\n");
	Free();
	delete WaveC;
	delete LevelC;
	delete PotC;
	delete RotPotC;
	delete[] Pot;
	//printf("Ende ~WaveFuncPlot\n");
}

void WaveFuncPlot::Free()
{
	if (WF != 0)
	{
		Destroy(WF, NChannels, Nv);
		delete[] E;
		if (P != 0)
		{
			delete[] P;
			delete[] N;
			delete[] F;
			delete[] ri;
			delete[] ra;
			P=0;
		}
		WF = 0;
	}
}

void WaveFuncPlot::PotChanged(QString PName)
{
	//printf("WavefuncPlot::PotChanged\n");
	if (MW == 0) return;
	Free();
	Pot[0] = MW->getPotential(PName);
	if (Pot[0] == 0) return;
	Molecule *mol = Pot[0]->getMolecule();
	if (mol == 0) return;
	int i, n, iso, N, *C;
	QString IsoN;
	if ((SWF = Pot[0]->getCoupledSineWaveFunc()) != 0)
	{
		//ElState *St;
		nDatenS = NChannels = SWF->getNChannels();
		Pot[0]->getCoupledPotentials(i, Pot, C);
	}
	else nDatenS = NChannels = 1;
	IsoTab *Iso = mol->getIso();
	PMax = Pot[0]->getPoint(Pot[0]->getMaxR());
	PMin = Pot[0]->getMinimumE();
	if (Mol == mol) 
	{
		if ((N = IsoBox->count()) == Iso->numIso) iso = IsoBox->currentIndex();
		else 
		{
			for (n=0, IsoN = IsoBox->currentText(); (n < N ? IsoN != Iso->getIsoName(n) : false); n++) ;
			iso = (n<N ? n : Iso->refIso);
		}
		if (SWF == 0)
		{
			if (IsoBox->currentText().toInt() == iso || N < Iso->numIso)
			{
				IsoBox->clear();
				for (i=0; i < Iso->numIso; i++) IsoBox->addItem(Iso->getIsoName(i));
				IsoBox->setCurrentIndex(iso);
			}
			JBox->setEditable(true);
		}
		else
		{
			int *I;
			SWF->getIso(N, I);
			IsoBox->clear();
			for (n=0; n<N; n++)
			{
				IsoBox->addItem(Iso->getIsoName(I[n]));
				if (I[n] == iso) IsoBox->setCurrentIndex(n);
			}
			JBox->blockSignals(true);
			JBox->setEditable(false);
			JBox->blockSignals(false);
		}
	}
	else
	{
		if (SWF == 0)
		{
			Mol = mol;
			IsoBox->clear();
			for (i=0; i < Iso->numIso; i++) IsoBox->addItem(Iso->getIsoName(i));
			IsoBox->setCurrentIndex(Iso->refIso);
			JBox->setEditable(true);
		}
		else
		{
			int *I;
			SWF->getIso(N, I);
			IsoBox->clear();
			for (n=0; n<N; n++)
			{
				IsoBox->addItem(Iso->getIsoName(I[n]));
				if (I[n] == Iso->refIso) IsoBox->setCurrentIndex(n);
			}
			JBox->blockSignals(true);
			JBox->setEditable(false);
			JBox->blockSignals(false);
		}
	}
	//WaveFuncSF = 1e-3;
	IsoChanged(IsoBox->currentIndex());
	//printf("WaveFuncSF=%g, E[0]=%g, PMin=%g, WMax=%g\n", WaveFuncSF, E[0], PMin, WMax);
	XMin = Pot[0]->getMinR();
	XMax = Pot[0]->getMaxR();
	YMin = PMin - 0.1 * (PMax - PMin);
	YMax = 2.0 * PMax - PMin;
	xStart->setText(QString::number(XMin));
    xStop->setText(QString::number(XMax));
    yStart->setText(QString::number(YMin, 'g', 5));
    yStop->setText(QString::number(YMax, 'g', 5));
	//Paint();
	//printf("Ende PotChanged\n");
}

void WaveFuncPlot::IsoChanged(int I)
{
	if (SWF != 0)
	{
		int J, NJ = SWF->getNJ(), cJ = JBox->currentText().toInt();
		if (cJ < 0) cJ = 0;
		JBox->blockSignals(true);
		JBox->clear();
		for (J=0; J < NJ; J++) if (SWF->isJA(I, J)) 
		{
			JBox->addItem(QString::number(J));
			if (cJ >= 0 && J >= cJ) 
			{
				JBox->setCurrentIndex(JBox->count() - 1);
				cJ = -1;
			}
		}
		JBox->blockSignals(false);
	}
	else 
	{
		if (JBox->count() < 336)
		{
			int cJ = JBox->currentText().toInt(), n;
			JBox->blockSignals(true);
			JBox->clear();
			for (n=0; n <= 335; n++) JBox->addItem(QString::number(n));
			if (cJ < 0) cJ = 0;
			if (cJ <= 335) JBox->setCurrentIndex(cJ);
			else
			{
				JBox->addItem(QString::number(cJ));
				JBox->setCurrentIndex(336);
			}
			JBox->blockSignals(false);
		}
	}
	JChanged();	
}

void WaveFuncPlot::JChanged()
{
	if (MW == 0) return;
	if (Pot[0] == 0) return;
	QString JText = JBox->currentText();
	if (JText.isEmpty()) return;
	int i, n, v = vBox->currentText().toInt(), J = JText.toInt(), Iso = IsoBox->currentIndex();
	if (Iso < 0) Iso = 0;
	Free();
	vBox->blockSignals(true);
	vBox->clear();
	if (SWF == 0)
	{
		WF = new double**[1];
        Pot[0]->getWaveFunc(Iso, J, Nv = cMaxv, E, WF[0], ri, ra, P, F, N, NumWFPoints);
		if (Nv == 0)
		{
			printf("No vibrational levels found!\n");
			return;
		}
		WaveFuncSF = (PMax - PMin) / (5.0 * WF[0][0][P[0]] * N[0]);
		//printf("P[0]=%d\n", P[0]);
		for (i=0; i < Nv; i++) vBox->addItem(QString::number(i));
		if (v >= Nv) v = Nv - 1;
		vBox->setCurrentIndex(v);
	}
	else
	{
		int *vs, j;
		double M = 0.0, m;
		if (J<0 && J >= SWF->getNJ()) return;
		SWF->getv(Iso, J, n, vs);
		if (n==0) return;
		if (v >= (Nv = vs[n - 1] + 1)) v = Nv - 1;
		for (i=0; i<n; i++)
		{
			vBox->addItem(QString::number(vs[i]));
			if (vs[i] == v) vBox->setCurrentIndex(i);
		}
		WF = Create(NChannels, Nv, NumWFPoints);
		E = new double[Nv];
		SWF->getWaveFuncs(Iso, J, rmin, rmax, NumWFPoints, WF, E);
		for (i=0, j = vs[0]; i < NChannels; i++) for (n=0; n < NumWFPoints; n++) if ((m = fabs(WF[i][j][n])) > M) M = m; 
		WaveFuncSF = (PMax - PMin) / (5.0 * M);
		delete[] vs;
	}
	vBox->blockSignals(false);
	vChanged(v);
}

void WaveFuncPlot::vChanged(int v)
{
	//printf("WaveFuncPlot::vChanged(%d), Nv=%d\n", v, Nv);
	int n;
	for (n=0; n < nDatenS; n++) Daten[n].reinit();
	if (WF != 0 && Nv > 0) 
	{
		if (v==-1) v=0;
		int i;
		double Max = 0.0, Min = 0.0;
		double R, V;
		if (SWF == 0)
		{
			//printf("N[%d]=%g, WaveFuncSF=%g\n", v, N[v], WaveFuncSF);
			//printf("ri[v]=%d, ra[v]=%d, P[v]=%d, NumWFPoints=%d\n", ri[v], ra[v], P[v], NumWFPoints);
			Ra = rmin + h * double(ri[v]);
			Ro = rmin + h * double(ra[v] - 1);
			//printf("Ra=%f, Ro=%f\n", Ra, Ro);
			//printf("ri=%d, ra=%d, NumWFPoints=%d\n", ri[v], ra[v], NumWFPoints);
			for (i = ri[v], R = Ra; i < ra[v]; i++, R+=h) 
			{	
				V = E[v] + (i <= P[v] ? WaveFuncSF * WF[0][v][i] * N[v] : WaveFuncSF * F[v] * WF[0][v][i] * N[v]);
				/*if (i == P[v]) 
				printf("V=%f, R=%f, Scale=%g, WF[%d][%d]=%g, F[%d]=%f\n", 
					   V, R, WaveFuncSF, v, i, WF[v][i], v, F[v]);*/
				if (V > Max) Max = V;
				if (V < Min) Min = V;
				Daten->AddValue(R, V, false);
			}
		}
		else
		{
			Ra = rmin;
			Ro = rmax;
			for (i=0, R = Ra; i < NumWFPoints; i++, R+=h) for (n=0; n < NChannels; n++)
			{
				V = E[v] + WaveFuncSF * WF[n][v][i];
				if (V > Max) Max = V;
				if (V < Min) Min = V;
				Daten[n].AddValue(R, V, false);
			}
		}
		if (xStart->text().toDouble() > Ra) 
			xStart->setText(QString::number(Ra - 0.1 * (Ro - Ra), 'g', 5));
		if (xStop->text().toDouble() < Ro) xStop->setText(QString::number(Ro + 0.1 * (Ro - Ra), 'g', 5));
		if (yStart->text().toDouble() > Min) 
			yStart->setText(QString::number(Min - 0.1 * (Max - Min), 'g', 5));
		if (yStop->text().toDouble() < Max) 
			yStop->setText(QString::number(Max + 0.1 * (Max - Min), 'g', 5));
	}
	Paint();
}

void WaveFuncPlot::PSpektrum(QPainter &P, const QRect &A, bool PrintFN)
{
	//printf("WaveFuncPlot::PSpektrum\n");
	int i, s, w = A.width() - ScaleYWidth, l = A.left() + ScaleYWidth - 1;
	int t = A.top(), lp, J = JBox->currentText().toInt(), Iso = IsoBox->currentIndex(), v = vBox->currentIndex();
	int p=0, b = A.bottom() - ScaleXHeight, j, r = l + w, LE, n;
    double minR = xStart->text().toDouble(), maxR = xStop->text().toDouble();
	double rpot[w+1], *pot, H = (maxR - minR) / double(w+1);
	if (v==-1 && Nv > 0) 
	{
		v=0;
		vBox->setCurrentIndex(0);
	}
	for (n=0; n < NChannels; n++)
	{
		pot = Pot[n]->getPoints(minR, maxR, w+1);
		if (pot == 0) continue;
		for (s=0; (s<=w ? pot[s] == 0.0 : false); s++) ;
		Pot[n]->CRotPot(J, Iso, pot, rpot, minR + H * double(s), H, w+1-s);
		for (j=0; j<2; j++)
		{
			if (j==0) 
			{
				if (NChannels == 1) P.setPen(*PotC);
				else switch (n)
				{
					case 0:
						P.setPen(QColor(0, 0, 0));
						break;
					case 1:
						P.setPen(QColor(255, 0, 0));
						break;
					case 2:
						P.setPen(QColor(0, 255, 0));
						break;
					case 3:
						P.setPen(QColor(0, 0, 255));
						break;
				}
			}
			else
			{
				delete[] pot;
				pot = rpot;
				if (NChannels == 1) P.setPen(*RotPotC);
				else switch (n)
				{
					case 0:
						P.setPen(QColor(100, 100, 100));
						break;
					case 1:
						P.setPen(QColor(255, 127, 127));
						break;
					case 2:
						P.setPen(QColor(127, 255, 127));
						break;
					case 3:
						P.setPen(QColor(127, 127, 255));
						break;
				}
			}
			for (i=s+l; i<=l+w; i++)
			{
				lp = p;
				p = int(YO - YSF * pot[i-l]);
				if (p > b) p = b;
				if (p < t) p = t;
				if (!(p==t && lp==t) && !(p==b && lp==b))
				{
					if (i>l+s ? (p > lp ? p - lp : lp - p) > 1 : false) P.drawLine(i-1, lp, i, p);
					else P.drawPoint(i, p);
				}
			}
		}
	}
	//printf("Nach Paint pots\n");
	if (v==-1) return;
	LE = int(YO - YSF * E[v]);
	p = int(XO + Ra * XSF);
	lp = int(XO + Ro * XSF);
	if (LE > t && LE < b && p<r && lp > l)
	{
		P.setPen(*LevelC);
		P.drawLine((p>l ? p : l), LE, (lp < r ? lp : r), LE);
	}
	DiagWindow::PSpektrum(P, A, PrintFN);
}
