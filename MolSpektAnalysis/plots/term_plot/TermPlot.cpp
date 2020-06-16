//
// C++ Implementation: TermPlot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "TermPlot.h"
#include "MainWindow.h"
#include "elstate.h"
#include "molecule.h"
#include "termtable.h"
#include "linetable.h"
#include "Spektrum.h"
#include "utils.h"
#include "tools.h"
#include "fitdata.h"
#include "isotab.h"
#include "tpviewlist.h"
#include "tableline.h"

#include <qlabel.h>
#include <qpushbutton.h>
#include <qfile.h>
#include <QTextStream>
#include <qstring.h>
#include <QFileDialog>
#include <qpainter.h>
#include <qcolor.h>
#include <qwidget.h>
#include <QGridLayout>

#include <math.h>

TermPlot::TermPlot(MainWindow *MW, Molecule *Mol, ElState *State) : DiagWindow(TermEnergyPlot, MW)
{
    //printf("Termplot::Termplot\n");
    int i;
	vs = NULL;
	iso = 0;
	mol = Mol;
	state = State;
	lstate = 0;
	LLT = 0;
	vss = Jss = 0;
	if (Mol != 0)
	{
		if (state != 0) 
			setWindowTitle("TermPlot of state " + state->getName() + " of the molecule " 
					+ Mol->getName());
		else setWindowTitle("TermPlot of an unknown state of the molecule " + Mol->getName());
	}
	else setWindowTitle("Unassigned termPlot");
	XUnit = "J(J+1)-\\Omega^2";
	YUnit = "Energy [cm^{-1}]";
	AnzahlLinien = AnzahlMarked = 0;
	LTE = MTE = 0;
	LJtJp1 = MJtJp1 = 0;
	MIso = LIso = 0;
	//printf("TermPlot: Nach ersten Initialisierungen\n");
	SMarked = false;
	LMarked = 0;
	if (mol != 0) 
	{
		IsoTab *IsoT = mol->getIso();
		for (i=0; i<IsoT->numIso; i++)
			IsoBox->addItem(QString::number(IsoT->mNumIso1[i]) + *IsoT->chSymb1 
					+ QString::number(IsoT->mNumIso2[i]) + *IsoT->chSymb2);

	}
	
    xStartLabel->setText("J*(J+1):    von ");
    xStart->setText("0");
    xStopLabel->setText(" bis ");
    xStop->setText("17000");
    yStartLabel->setText("Energie in cm-1:    von ");
    yStart->setText("18000");
    yStopLabel->setText(" bis ");
    yStop->setText("21000");
    //Durchsuchen->setGeometry(680, 10, 100, 20);
    setMinimumSize(790, 600);
    connect(Evs, SIGNAL(returnPressed()), this, SLOT(Updatevs()));
	connect(Bild, SIGNAL(ScreenSplitted(int, int, int, int)), 
			this, SLOT(SMark(int , int, int, int)));
	connect(IsoBox, SIGNAL(currentIndexChanged(int)), this, SLOT(IsoChanged(int)));
    //printf("TermPlot: Vor Update()\n");
    UpdateLines();
	UpdateMarked();
}

TermPlot::~TermPlot()
{
	if (vs != NULL) delete[] vs;
    if (LTE != NULL) delete[] LTE;
	if (MTE != 0) delete[] MTE;
    if (LJtJp1 != NULL) delete[] LJtJp1;
	if (MJtJp1 != 0) delete[] MJtJp1;
	if (LMarked != NULL) delete[] LMarked;
	if (MIso != 0) delete[] MIso;
	if (LIso != 0) delete[] LIso;
	if (LLT != 0) delete[] LLT;
	if (vss != 0) delete[] vss;
	if (Jss != 0) delete[] Jss;
}

/*void TermPlot::Mark()
{
	bool M = !Bild->splitScreenEnabled();
	Bild->setSplitScreen(M);
	editMenu->setItemChecked(2, M);
	if (M) editMenu->setItemChecked(1, false);
	else SMarked = false;
}*/

void TermPlot::Zoom()
{
	bool Z = !Bild->ZoomEnabled();
	Bild->setZoom(Z);
	//editMenu->setItemChecked(1, Z);
	//if (Z) editMenu->setItemChecked(2, false);
}

/*void TermPlot::SMark(int x1, int y1, int x2, int y2)
{
	printf("SMark: x1=%d, y1=%d, x2=%d, y2=%d\n", x1, y1, x2, y2);
	slx1 = x1;
	sly1 = y1;
	slx2 = x2;
	sly2 = y2;
	SMarked = true;
	Paint();
}

void TermPlot::PlotDunham(int KMax, int *LMax, double **Par, double **Korr)
{
	kMax = KMax;
    lMax = LMax;
    DunhamP = Par; 
    DKorr = Korr;
	Paint();
}*/

void TermPlot::Updatevs()
{
    /*int i, v = Evs->text().toInt();
    for (i=0; i<AnzahlMarker; i++) 
		if (marker[i].Marked && marker[i].DisplayData) marker[i].vs = v;
	if (SMarked) 
	{
		for (i=0; i<AnzahlLinien; i++) if (LMarked[i]) vs[i] = v;
		(**tabelle).Updatevs(vs);
	}
    emit vsUpdated();*/
}

void TermPlot::IsoChanged(int nIso)
{
	iso = nIso;
	Paint();
}

void TermPlot::SetBoundaries()
{
	//printf("Termplot::SetBoundaries()\n");
	int i;
    XMin = XMax = 0.0;
    YMin=50000.0, YMax = 0.0;
    if (state != 0 ? state->getTermTable() != 0 : false)
	{
		int Omega = state->getLambda();
		int MJ = state->getTermTable()->getMaxJ();
		XMax = MJ * MJ + 1 - Omega * Omega;
		double ****D = state->getTermTable()->getData();
		int NC = state->getTermTable()->getNumComp();
		int MV = state->getTermTable()->getMaxv();
		int NI = state->getTermTable()->getNumIso();
		int v, I, C, J;
		for (C=0; C<NC; C++) for (I=0; I<NI; I++) for (v=0; v<=MV; v++) for (J=0; J<=MJ; J++)
		{
			if (D[C][I][v][J] < YMin) YMin = D[C][I][v][J];
			if (D[C][I][v][J] > YMax) YMax = D[C][I][v][J];
		}
	}
	//printf("XMax=%g, YMin=%g\n", XMax, YMin);
	for (i=0; i < AnzahlMarked; i++) if (MTE[i] != -1.0)
    {
		if (MJtJp1[i] > XMax) XMax = MJtJp1[i];
		if (MTE[i] < YMin) YMin = MTE[i];
		else if (MTE[i] > YMax) YMax = MTE[i];
    }
	//printf("XMax=%g, YMin=%g\n", XMax, YMin);
    for (i=0; i<AnzahlLinien; i++) if (LTE[i] != -1.0)
    {
		//if (MTE[i] < 1e4) printf("LTE[%d]=%g\n", i, LTE[i]);
		if (LJtJp1[i] > XMax) XMax = LJtJp1[i];
		if (LTE[i] < YMin) YMin = LTE[i];
		else if (LTE[i] > YMax) YMax = LTE[i];
    }
	//printf("XMax=%g, YMin=%g\n", XMax, YMin);
    XMin -= 0.05 * (XMax - XMin);
	XMax += 20.0;
    YMin -= 50.0;
    YMax += 50.5;
    xStart->setText(QString::number(XMin));
    xStop->setText(QString::number(XMax));
    yStart->setText(QString::number(int(YMin)));
    yStop->setText(QString::number(int(YMax)));
	//printf("XMax=%g, YMin=%g\n", XMax, YMin);
    Paint();
	//printf("Ende setBoundaries\n");
}

void TermPlot::drawP(int x, int y, QColor C, QPainter &p)
{
    p.setPen(C);
    p.drawLine(x - 3, y, x + 3, y);
    p.drawLine(x, y - 3, x, y + 3);
}

void TermPlot::PSpektrum(QPainter &p, const QRect & A, bool /*PrintFN*/ )
{
    //printf("Termplot::PSpektrum()\n");
    int i, j=0, J, l = ScaleYWidth, t = 1, NLTSel, *LTSelJs, n; 
	int h = A.height() - t - ScaleXHeight + 1, w = A.width() - l + 1/*, aEp*/;
    int minJ = (int)xStart->text().toDouble(), maxJ = (int)xStop->text().toDouble()/*, v*/;
	//int LPos[AnzahlLinien][2];
    double minE = yStart->text().toDouble(), maxE = yStop->text().toDouble()/*, TermE[2000], vj, vf*/;
    //double Ef, Ee, Eft, Eet;
    int Omega = (state != 0 ? state->getOmega() : 0);
	if (maxE - minE < 1.0) yStop->setText(QString::number(maxE = minE + 1.0));
    if (maxJ - minJ < 1) xStop->setText(QString::number(maxJ = minJ + 1));
    double JSc = (double)w / (maxJ - minJ), ESc = (double)(h - t) / (maxE - minE), *LTSelE;
    //printf("j=%d, maxJ=%d, minE=%f, maxE=%f, JSc=%f, ESc=%f\n", j, maxJ, minE, maxE, JSc, ESc);
    int Jp, Ep, NFData;
	Transition *T;
	LineTable *L;
	FitData *FD;
	TableLine *FData;
	TermTable *TT = 0;
	QList<TPViewList> SList, VList;
	TPViewList CSel, CView;
	if (state != 0 ? (TT = state->getTermTable()) != 0 : false)
	{
		int MJ = TT->getMaxJ();
		XMax = MJ * MJ + 1 - Omega * Omega;
		double ****D = TT->getData();
		int NC = TT->getNumComp();
		int MV = TT->getMaxv();
		int v, C, J;
		for (J=0; J<=MJ; J++) 
		{
			Jp = l + int((J*(J+1) - Omega * Omega - minJ) * JSc);
			for (C=0; C<NC; C++) for (v=0; v<=MV; v++)
			{
				if (D[C][iso][v][J] == 0.0 || (v>0 ? D[C][iso][v][J] < D[C][iso][v-1][J] : false)
				    || (J>0 ? D[C][iso][v][J] < D[C][iso][v][J-1] : false)) break; 
				Ep = t + int((maxE - D[C][iso][v][J]) * ESc);
				if (Jp >= l && Ep >= t && Jp <= l+w && Ep <= t+h) 
					drawP(Jp, Ep, QColor(175, 175, 175), p);
			}
		}
	}
	//printf("Vor Linien\n");
    for (i=0; i<AnzahlLinien; i++) if (LTE[i] != -1.0 && LIso[i] == iso)
	{
		Jp = /*LPos[i][0] =*/ l + (int)((LJtJp1[i] - Omega * Omega - minJ) * JSc);
		Ep = /*LPos[i][1] =*/ t + (int)((maxE - LTE[i]) * ESc);
		if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(0, 0, 200), p);
		//printf("i=%d, Jp=%d, Ep=%d\n", i, Jp, Ep);
	} //else printf("LTE[%d]=%g, LIso[%d]=%d, iso=%d\n", i, LTE[i], i, LIso[i], iso);}
	//printf("Vor marked\n");
	if (state != 0) for (n=0; n < state->getNumFitDataSets(); n++) if ((FD = state->getFitData(n)) != 0)
	{
		FD->getData(FData, NFData);
		for (i = CSel.NP = CView.NP = 0; i < NFData; i++) 
			if (FData[i].Iso == iso && FData[i].isTE)
		{
			if (FData[i].isSelected) CSel.NP++;
			if (FData[i].isViewn) CView.NP++;
		}
		if (CSel.NP > 0)
		{
			CSel.EP = new int[CSel.NP];
			CSel.JP = new int[CSel.NP];
		}
		if (CView.NP > 0)
		{
			CView.EP = new int[CView.NP];
			CView.JP = new int[CView.NP];
		}
		for (i = CSel.NP = CView.NP = 0; i < NFData; i++) 
			if (FData[i].Iso == iso && FData[i].isTE)
		{
			Jp = l + (int)((FData[i].Jss * (FData[i].Jss + 1) - Omega * Omega - minJ) * JSc);
			Ep = t + (int)((maxE - FData[i].WN) * ESc);
			if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h)
			{
				drawP(Jp, Ep, QColor(0, 0, 200), p);
				if (FData[i].isSelected) 
				{
					CSel.EP[CSel.NP] = Ep;
					CSel.JP[CSel.NP++] = Jp;
				}
				if (FData[i].isViewn)
				{
					CView.EP[CView.NP] = Ep;
					CView.JP[CView.NP++] = Jp;
				}
			}
		}
		if (CSel.NP > 0) SList.append(CSel);
		if (CView.NP > 0) VList.append(CView);
	}
	if (TT != 0)
	{
		TT->getViewnE(LTSelJs, LTSelE, NLTSel);
		if (NLTSel > 0)
		{
			for (i=0; i < NLTSel; i++)
			{
				Jp = l + (int)((LTSelJs[i] * (LTSelJs[i] + 1) - Omega * Omega - minJ) * JSc);
				Ep = t + (int)((maxE - LTSelE[i]) * ESc);
				if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(0, 200, 0), p);
			}
			delete[] LTSelE;
			delete[] LTSelJs;
		}
	}
	if (mol != 0) 
	{
		for (j=0; j < mol->getNumTransitions(); j++) if ((T = mol->getTransitionP(j))->getUpperState() == state)
			for (J=0; J < T->getNumLineTables(); J++)
		{
			L = T->getLineTable(J);
			if (L == 0) continue;
			L->getViewnE(LTSelJs, LTSelE, NLTSel);
			if (NLTSel > 0)
			{
				for (i=0; i < NLTSel; i++)
				{
					Jp = l + (int)((LTSelJs[i] * (LTSelJs[i] + 1) - Omega * Omega - minJ) * JSc);
					Ep = t + (int)((maxE - LTSelE[i]) * ESc);
					if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(0, 200, 0), p);
				}
				delete[] LTSelE;
				delete[] LTSelJs;
			}
		}
	}
	for (n=0; n < VList.count(); n++) for (i=0; i < VList[n].NP; i++) drawP(VList[n].JP[i], VList[n].EP[i], QColor(0, 200, 0), p);
	if (TT != 0)
	{
		TT->getSelE(LTSelJs, LTSelE, NLTSel);
		if (NLTSel > 0)
		{
			for (i=0; i < NLTSel; i++)
			{
				Jp = l + (int)((LTSelJs[i] * (LTSelJs[i] + 1) - Omega * Omega - minJ) * JSc);
				Ep = t + (int)((maxE - LTSelE[i]) * ESc);
				if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(150, 150, 0), p);
			}
			delete[] LTSelE;
			delete[] LTSelJs;
		}
	}
	if (mol != 0) 
	{
		for (j=0; j < mol->getNumTransitions(); j++) if ((T = mol->getTransitionP(j))->getUpperState() == state)
			for (J=0; J < T->getNumLineTables(); J++)
		{
			L = T->getLineTable(J);
			if (L == 0) continue;
			L->getSelData(LTSelJs, LTSelE, NLTSel);
			for (i=0; i < NLTSel; i++)
			{
				Jp = l + (int)((LTSelJs[i] * (LTSelJs[i] + 1) - Omega * Omega - minJ) * JSc);
				Ep = t + (int)((maxE - LTSelE[i]) * ESc);
				if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(150, 150, 0), p);
			}
		}
	}
	else 
		QMessageBox::information(this, "MolSpektAnalysis", 
						"This plot might have contained more information if the table had been assigned to a molecule first.");
	for (n=0; n < SList.count(); n++) for (i=0; i < SList[n].NP; i++) 
		drawP(SList[n].JP[i], SList[n].EP[i], QColor(150, 150, 0), p);
	for (i=0; i < AnzahlMarked; i++) if (MTE[i] != -1.0 && MIso[i] == iso)
	{
		Jp = l + (int)((MJtJp1[i] - Omega * Omega - minJ) * JSc);
		Ep = t + (int)((maxE - MTE[i]) * ESc);
		if (Jp >= l && Ep >= t && Jp <= l + w && Ep <= t + h) drawP(Jp, Ep, QColor(200, 0, 0), p);
	}
	/*if (SMarked) 
	{
		double St = (double)(sly2 - sly1) / (slx2 - slx1);
		for (i=0; i<AnzahlLinien; i++) 
		{
			if ((LPos[i][0] >= slx2 || (LPos[i][0] >= slx1 
						 && LPos[i][1] >= sly1 + (int)(St * (LPos[i][0] - slx1)))) && vs[i] == -1)
			{
				drawP(LPos[i][0], LPos[i][1], QColor(0, 0, 255), p);
				LMarked[i] = true;
			}
			else LMarked[i] = false;
		}
	}
	
    if (lMax != NULL) for (v=0, vf=0.5; Jp > 1 || v==0; v++)
    {
		printf("Beginn Zeichnen der Dunhamkoeffizienten, v=%d\n", v);
		for (Jp=0, Ep = h; Ep >= 0 && Jp <= w && (Jp==0 || Ep <= aEp); Jp++)
		{
	    	vj = double(Jp) / JSc + minJ - omegaSq;
	    	Ee = Ef = 0.0;
	    	for (i=kMax; i>=0; i--) 
	    	{
				Eft = Eet = 0.0;
				for (j=lMax[i]; j>=0; j--)
				{
		    		Eet *= vj;
		    		Eft *= vj;
		    		Eet += DunhamP[i][j] + DKorr[i][j];
		    		Eft += DunhamP[i][j];
				}
				Ee *= vf;
				Ef *= vf;
				Ee += Eet;
				Ef += Eft;
	    	}
			aEp = Ep;
			if ((Ep = h - (int)((Ee - minE) * ESc)) < h)
	    	{
				if (v==0) printf("Ep=%d, h=%d, Jp=%d, w=%d\n", Ep, h, Jp, w);
				p.setPen(QColor(255, 0, 0));
	        	p.drawPoint(Jp, Ep);
	    	}
	    	if ((Ep = h - (int)((Ef - minE) * ESc)) < h)
	 		{
				if (v==0) printf("Ep=%d, Jp=%d\n", Ep, Jp);
	        	p.setPen(QColor(0, 255, 0));
	        	p.drawPoint(Jp, Ep);
	 		}
	 		if (v==0) printf("vj=%f, Jp=%d, Ep=%d, h=%d\n", vj, Jp, Ep, h);
		}
		vf += 1.0;
   	}*/
	//printf("Ende von PSpektrum\n");
}

void TermPlot::UpdateMarked()
{
	if (ELU == 0 || MW == 0) return;
	if (MTE != 0) delete[] MTE;
	if (MJtJp1 != 0) delete[] MJtJp1;
	if (MIso != 0) delete[] MIso;
	int n, m, S, NS = MW->getNumSpectra(), AM;
	Marker *marker;
	TermTable *TT;
	for (S=AnzahlMarked=0; S < NS; S++) 
	{
		MW->getSpectrum(S)->GetMarker(AM, marker);
		for (m=0; m < AM; m++)
			if (marker[m].Marked && marker[m].DisplayData && marker[m].UState == state) AnzahlMarked++;
	}
	if (AnzahlMarked == 0)
	{
		MTE = 0;
		MJtJp1 = 0;
		MIso = 0;
		return;
	}
	for (S=n=0; S < NS; S++)
	{
		MTE = new double[AnzahlMarked];
		MJtJp1 = new int[AnzahlMarked];
		MIso = new int[AnzahlMarked];
		MW->getSpectrum(S)->GetMarker(AM, marker);
		for (m=0; m < AM; m++)
			if (marker[m].Marked && marker[m].DisplayData && marker[m].UState == state)
		{
			if (marker[m].LState != lstate)
			{
				if ((TT = ((lstate = marker[m].LState) != 0 ? lstate->getTermTable() : 0)) != 0)
				{
					ELU = TT->getData();
					nlC = TT->getNumComp();
					nlI = TT->getNumIso();
					nlv = TT->getMaxv() + 1;
					nlJ = TT->getMaxJ() + 1;
				}
				else
				{
					ELU = 0;
					nlC = nlI = nlJ = nlv = 0;
				}
			}
			if (marker[m].vss >= 0 && marker[m].vss < nlv
					 && marker[m].Iso >= 0 && marker[m].Iso < nlI && marker[m].Jss >= 0
					 && marker[m].Jss < nlJ)
			{
				MIso[n] = marker[m].Iso;
				MTE[n] = marker[m].Line[0] + ELU[0][marker[m].Iso][marker[m].vss][marker[m].Jss];
				MJtJp1[n++] = marker[m].Js * (marker[m].Js + 1);
			}
		}
	}
	AnzahlMarked = n;
	SetBoundaries();
}

void TermPlot::UpdateLines()
{
    //printf("Termplot::UpdateLines()\n");
    if (LTE != NULL) delete[] LTE;
    if (LJtJp1 != NULL) delete[] LJtJp1;
	if (vs != NULL) delete[] vs;
	if (LMarked != NULL) delete[] LMarked;
	if (LIso != 0) delete[] LIso;
	if (LLT != 0) delete[] LLT;
	if (vss != 0) delete[] vss;
	if (Jss != 0) delete[] Jss;
	//printf("Nach dem Löschen\n");
	int i, j, n, l, li, numLineTables, numTransitions = (mol != 0 ? mol->getNumTransitions() : 0);
	Transition *T;
	LineTable *L;
	TermTable *TT;
	for (n = AnzahlLinien = 0; n < numTransitions; n++) 
		if ((T = mol->getTransitionP(n))->getUpperState() == state)
			for (numLineTables = T->getNumLineTables(), l=0; l < numLineTables; l++) if ((L = T->getLineTable(l)) != 0)
				AnzahlLinien += L->getAnzahlLinien();
	//printf("Vor Hauptteil\n");
	if (AnzahlLinien > 0)
	{
		LTE = new double[AnzahlLinien];
		LJtJp1 = new int[AnzahlLinien];
		vs = new int[AnzahlLinien];
		LMarked = new bool[AnzahlLinien];
		LIso = new int[AnzahlLinien];
		vss = new int[AnzahlLinien];
		Jss = new int[AnzahlLinien];
		LLT = new LineTable*[AnzahlLinien];
		double *Unc = new double[AnzahlLinien];
		int **LZuordnung = CreateInt(AnzahlLinien, 6);
		//printf("Nach Create\n");
		for (i = n = 0; n < numTransitions; n++) 
			if ((T = mol->getTransitionP(n))->getUpperState() == state)
		{
			if (T->getLowerState() != lstate)
			{
				lstate = T->getLowerState();
				if ((TT = lstate->getTermTable()) != 0) 
				{
					ELU = TT->getData();
					nlC = TT->getNumComp();
					nlI = TT->getNumIso();
					nlv = TT->getMaxv() + 1;
					nlJ = TT->getMaxJ() + 1;
				}
				else
				{
					ELU = 0;
					nlC = nlI = nlv = nlJ = 0;
				}
			}
			for (numLineTables = T->getNumLineTables(), l=0; l < numLineTables; l++)
			{
				//printf("i=%d, n=%d\n", i, n);
				li = i;
				L = T->getLineTable(l);
				if (L==0) continue;
				L->getLines(LZuordnung + i, LTE + i, Unc + i);
				i += L->getAnzahlLinien();
				for (j = li; j < i; j++) 
				{
					if (LZuordnung[j][0] >= 0 && LZuordnung[j][0] < nlI && LZuordnung[j][3] >= 0 
							&& LZuordnung[j][3] < nlv && LZuordnung[j][4] >= 0 && LZuordnung[j][4] < nlJ)
						LTE[j] += ELU[0][LZuordnung[j][0]][LZuordnung[j][3]][LZuordnung[j][4]];
					else 
					{
						LTE[j] = -1.0;
						//printf("LZuordnun[%d][0]=%d, nlI=%d, LZuordnung[%d][3]=%d, nlJ=%d, LZuordnung[%d][4]=%d, nlv=%d\n", j, LZuordnung[j][0], nlI, j, LZuordnung[j][3], nlJ, j, LZuordnung[j][4], nlv);
					}
					LIso[j] = LZuordnung[j][0];
					LJtJp1[j] = LZuordnung[j][2] * (LZuordnung[j][2] + 1);
					vs[j] = LZuordnung[j][1];
					vss[j] = LZuordnung[j][3];
					Jss[j] = LZuordnung[j][4];
					LLT[j] = L;
					//printf("LTE[%d]=%f, LJtJp1[%d]=%f\n", i, LTE[i], i, LJtJp1[i]);
				}
			}
		}
		//printf("Vor Destroy\n");
		Destroy(LZuordnung, AnzahlLinien);
		delete[] Unc;
	}
	else
	{
		LTE = 0;
		LJtJp1 = 0;
		vs = 0;
		LMarked = 0;
		LIso = 0;
		vss = 0;
		Jss = 0;
		LLT = 0;
    }
    //printf("Ende von UpdateLines\n");
    SetBoundaries();
}

void TermPlot::PictureClicked(QPoint* P)
{
    if (AnzahlLinien == 0 || mol == 0) return;
	int minJ = (int)xStart->text().toDouble(), maxJ = (int)xStop->text().toDouble();
    double minE = yStart->text().toDouble(), maxE = yStop->text().toDouble();
    int Omega = (state != 0 ? state->getOmega() : 0);
	QRect A = Bild->contentsRect();
	double JSc = double(A.width() - ScaleYWidth + 1) / (maxJ - minJ);
	double ESc = (double)(A.height() - ScaleXHeight - 1) / (maxE - minE);
	int n, m, s=0, N=0, Jm = int(double(P->x() - ScaleYWidth - 3) / JSc) + Omega * Omega + minJ;
	int JM = int(double(P->x() - ScaleYWidth + 3) / JSc) + Omega * Omega + minJ, NI = 0, MaxJ = 0, Maxv = 0;
	double Em = maxE - double(P->y() + 2) / ESc, EM = maxE - double(P->y() - 4) / ESc;
	for (n=0; n < AnzahlLinien; n++) 
		if (LTE[n] != -1.0 && LTE[n] >= Em && LTE[n] <= EM && LJtJp1[n] >= Jm && LJtJp1[n] <= JM) N++;
	if (N==0) return;
	int I[N], v[N], J[N];
	double E[N], ****ELU = 0;
	LineTable *L=0;
	TermTable *T;
	Transition *Tr;
	ElState *lowState;
	for (n=m=0; n < AnzahlLinien; n++) 
		if (LTE[n] != -1.0 && LTE[n] >= Em && LTE[n] <= EM && LJtJp1[n] >= Jm && LJtJp1[n] <= JM)
	{
		if (LLT[n] != L)
		{
			if ((Tr = LLT[n]->getTransition()) != 0 ? ((lowState = Tr->getLowerState()) != 0 ?
				(T = lowState->getTermTable()) != 0 : false) : false)
			{
				ELU = T->getData();
				NI = T->getNumIso();
				MaxJ = T->getMaxJ();
				Maxv = T->getMaxv();
			}
			else NI = 0;
		}
		I[m] = LIso[n];
		v[m] = vs[n];
		J[m] = int(sqrt(double(LJtJp1[n])));
		E[m++] = (LIso[n] < NI && Jss[n] <= MaxJ && Maxv >= vss[n] ? LTE[n] - ELU[0][LIso[n]][vss[n]][Jss[n]] : 0.0);
		if (LLT[n] != L)
		{
			if (m>1) L->MarkLines(I+s, v+s, J+s, E+s, m-s-1);
			s=m-1;
			L = LLT[n];
		}
		if (m==N) L->MarkLines(I+s, v+s, J+s, E+s, m-s);
	}
}
