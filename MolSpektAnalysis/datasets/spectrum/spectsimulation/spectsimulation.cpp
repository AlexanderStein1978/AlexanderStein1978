//
// C++ Implementation: spectsimulation
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "spectsimulation.h"
#include "constants.h"
#include "MainWindow.h"
#include "molecule.h"
#include "elstate.h"
#include "utils.h"
#include "tools.h"
#include "potential.h"
#include "termtable.h"
#include "simline.h"
#include "isotab.h"
#include "doublevector.h"
#include "datensatz.h"

#include <math.h>

#include <QPainter>
#include <QComboBox>
#include <QLineEdit>
#include <QPaintDevice>
#include <QRect>
#include <QString>
#include <QTextStream>
#include <QStringList>
#include <QChar>
#include <QFile>
#include <QMessageBox>

SpectSimulation::SpectSimulation(MainWindow *mw) 
	: DiagWindow(SpectrumSimulation, mw, "Simulated spectra (*.sspk)", ".sspk")
{
	int i, N = (mw != 0 ? mw->getNumMolecules() : 0);
	//printf("SpektSimulation::SpektSimulation: NumMolecules=%d, mw=%d, MW=%d\n", N, mw, MW);
	Lines = 0;
	NumLines = 0;
	Iso = 0;
	Sim = 0;
	xStartLabel->setText("Frequency [cm-1]:    from ");
    xStart->setText("10000");
    xStopLabel->setText(" to ");
    xStop->setText("20000");
    yStartLabel->setText("Intensity:    from ");
    yStart->setText("0");
    yStopLabel->setText(" to ");
    yStop->setText("1000");
	SimBox->addItem("Absorption");
	SimBox->addItem("LIF");
	SimBox->addItem("Filtered absorption");
	for (i=0; i<N; i++) MolBox->addItem(MW->getMolecule(i)->getName());
	if (MW != 0) Mol = MW->getMolecule(0);
	TUnit->addItem(QString("°C").right(2));
	TUnit->addItem("K");
	TUnit->setEditable(false);
	LaserF->setEnabled(false);
	WindowE->setEnabled(false);
	WindowS->setEnabled(false);
	connect(Calc, SIGNAL(clicked()), this, SLOT(Simulate()));
	connect(MolBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(MolBoxChanged(QString)));
	connect(SimBox, SIGNAL(currentIndexChanged(int)), this, SLOT(SimBoxChanged(int)));
	printf("Ende SpektSimulation\n");
}

SpectSimulation::~SpectSimulation()
{
	if (Lines != 0) delete[] Lines;
	if (Iso != 0) delete Iso;
}

void SpectSimulation::cExcProfile(double Cf, double I, double Lf, double &Rs, double &Re, 
								  double *&Data, double RRes, double IRes, double DF, double L_gsq)
		 							const
{
	DoubleVector D;
	double RS = Lf - Cf, V = IRes, wdD, wdL;
	int i, is = 1;
	for (wdD = 0.0, i=0, wdL = RS; V >= IRes; i += is, wdD += RRes, wdL -= RRes)
	{
		D[i] = V = I * exp(DF * wdD * wdD) / sqrt(wdL * wdL + L_gsq);
		if (V < IRes && is == 1)
		{
			Re = Cf + wdD;
			i = is = -1;
			wdL = -RRes - RS;
			V = IRes;
		}
	}
	Rs = Cf -= wdD;
	Data = D.getData();
}

void SpectSimulation::flLine(double Cf, double Lf, double I, double &Rs, double &Re, double *&Data,
							 double RRes, double IRes, double L_gsq, double EPRs, double EPRe, 
		                     double *EPData) const
{
	int i, j, js=0;
	double DF = Cf / Lf;
	double RS, RSD = RRes * DF, V, wd, CfL, Ef = (EPRe - Lf) * DF + Cf;
	DoubleVector D;
	for (i=0, CfL = (EPRs - Lf) * DF + Cf, RS = Rs; CfL < Ef; i++, CfL += RSD)
	{
		while (RS < CfL) 
		{
			RS += RRes;
			js++;
		}
		for (j = js, wd = RS - CfL, V = IRes; V >= IRes; j++, wd += RRes) 
			D[j] += V = EPData[i] * I / sqrt(wd * wd + L_gsq);
		for (j = js - 1, wd = RS - CfL - RRes, V = IRes; V >= IRes; j--, wd -= RRes)
			D[j] += V = EPData[i] * I / sqrt(wd * wd + L_gsq);
	}
	Data = D.getData();
	Re = Rs + D.getLastIndex() * RRes;
	Rs += D.getFirstIndex() * RRes;
}

void SpectSimulation::getStates(Molecule *Mol, int &N, ElState **&ES) const
{
	int n, m, s;
	m = N = Mol->getNumStates();
	for (n=0; n<m; n++) if (Mol->getStateP(n)->getPotential() == 0) N--;
	ES = new ElState*[N];
	for (n=s=0; n<m; n++) if ((ES[s] = Mol->getStateP(n))->getPotential() != 0) s++;
}

void SpectSimulation::heapSortLines(SimLine *&Lines, int N) const
{
	SimLine B, *RB = new SimLine[N];
	int n, m, i;
	while (B.Sl != 0)
	{
		B.Sl = 0;
		for (n=0, m=2; m<=N; n++, m+=2) 
		{
			if (Lines[m-1].R < Lines[n].R)
			{
				B = Lines[m-1];
				Lines[m-1] = Lines[n];
				Lines[n] = B;
			}
			if (m < N ? Lines[m].R < Lines[n].R : false)
			{
				B = Lines[m];
				Lines[m] = Lines[n];
				Lines[n] = B;
			}
		}
	}
	for (i=0; i<N; i++)
	{
		RB[i] = Lines[0];
		for (n=0, m=2; Lines[n].Sl != 0; m=2*((n=m)+1))
		{
			if (m < N)
			{
				if (Lines[m-1].R < Lines[m].R) m--;
			}
			else if (m==N) m--;
			else break;
			Lines[n] = Lines[m];
		}
		if (m >= N) Lines[n].Sl = 0;
	}
	delete[] Lines;
	Lines = RB;
}

void SpectSimulation::MolBoxChanged(QString NMol)
{
	if (Mol != 0 ? NMol != Mol->getName() : true)
	{
		Mol = MW->getMolecule(NMol);
		if (Iso != 0) delete Iso;
		if (Mol != 0) Iso = Mol->getIso();
		else Iso = 0;
	}
}

void SpectSimulation::PSpektrum(QPainter &P, const QRect &R, bool PrintFN)
{
	printf("SpektSimulation::PSpektrum\n");
	DiagWindow::PSpektrum(P, R, PrintFN);
	if (NumLines == 0) return;
	int n, x, y1, y2 = (YO < YMax ? YO : YMax), th = TextHeight(QFont(), Lines[0].Su->getName());
	double XStart = xStart->text().toDouble(), XStop = xStop->text().toDouble();
	for (n=0; (n < NumLines ? Lines[n].R <= XStart : false); n++) ;
	for (P.setPen(QColor(0, 0, 255)); (n < NumLines ? Lines[n].R < XStop : false); n++) 
		if (Lines[n].I > XStart)
	{
		y1 = int(YO - YSF * (Lines[n].I));
		if (y1 < YMin) y1 = YMin;
		x = int(XO + XSF * Lines[n].R);
		P.drawLine(x, y2, x, y1);
		if ((n > 0 ? Lines[n].R - Lines[n-1]. R > th || Lines[n].I > Lines[n-1].I 
				   : Lines[n].R > XStart + th / 2 + 2)
		   && (n < NumLines - 1 ? Lines[n+1].R - Lines[n].R > th || Lines[n].I > Lines[n+1].I
			    				: Lines[n].R < XStop - th / 2 - 2)
		   && (y1 -= 20) > YMin)
			WriteText(P, x, y1, Iso->texName[Lines[n].Iso] + " " + Lines[n].Su->getName() 
					+ (Sim == 1 ? " \\rightarrow " : " \\leftarrow ") + QString("v'=") 
					+ QString::number(Lines[n].vs) + " J'=" + QString::number(Lines[n].Js) + " v''="
					+ QString::number(Lines[n].vss) + " J''=" + QString::number(Lines[n].Jss), 
									  QFont(), 1);
	}
	printf("Ende SpektSimulation::PSpektrum\n");
}

bool SpectSimulation::readData(QString FileName)
{
	int ST = 0, n, m, N, l=0, NI = Iso->numIso;
	QString Buffer, IsoN[NI];
	QStringList SList;
	QFile Datei(FileName);
	if (!read(&Datei)) return false;
	QTextStream S(&Datei);
	ElState *St;
	if ((Buffer = S.readLine()).left(19).toUpper() != "SIMULATED SPECTRUM:")
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), tr("Wrong file format!"));
		return false;
	}
	for (n=0; n < NI; n++) IsoN[n] = Iso->getIsoName(n);
	Daten->reinit();
	while (!S.atEnd())
	{
		Buffer = S.readLine();
		if (Buffer.indexOf("Table of lines", 0, Qt::CaseInsensitive) != -1)
		{
			ST = 1;
			continue;
		}
		if (Buffer.indexOf("Table of data points", 0, Qt::CaseInsensitive) != -1)
		{
			ST = 2;
			continue;
		}
		switch (ST)
		{
			case 0:
				if (Buffer.indexOf("Type", 0, Qt::CaseInsensitive) != -1)
				{
					if (Buffer.indexOf("Absorption", 0, Qt::CaseInsensitive) != -1)
						SimBox->setCurrentIndex(0);
					else if (Buffer.indexOf("LIF", 0, Qt::CaseInsensitive) != -1)
						SimBox->setCurrentIndex(1);
					else if (Buffer.indexOf("Filtered absorption", 0, Qt::CaseInsensitive) != -1)
						SimBox->setCurrentIndex(2);
					continue;
				}
				if ((n = Buffer.indexOf("Energy range from ", 0, Qt::CaseInsensitive)) != -1)
				{
					m = Buffer.indexOf("to", 0, Qt::CaseInsensitive);
					n+=18;
					xStart->setText(Buffer.mid(n, m-n-1));
					xStop->setText(Buffer.right(Buffer.length() - m - 3));
					continue;
				}
				if ((n = Buffer.indexOf("Intensity range form ", 0, Qt::CaseInsensitive)) != -1)
				{
					m = Buffer.indexOf("to", 0, Qt::CaseInsensitive);
					n+=21;
					yStart->setText(Buffer.mid(n, m-n-1));
					yStop->setText(Buffer.right(Buffer.length() - m - 3));
					continue;
				}
				if (Buffer.indexOf("Molecule", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":", 0) + 1; Buffer[n].isSpace(); n++) ;
					for (m = Buffer.length() -1; Buffer[m].isSpace(); m--) ;
					Buffer = Buffer.mid(n, m-n+1);
					for (n=0, N = MW->getNumMolecules(); 
										(n<N ? Buffer != MW->getMolecule(n)->getName() : false); n++) ;
					if (n<N) MolBox->setCurrentIndex(n);
					continue;
				}
				if (Buffer.indexOf("Temperature", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":", 0) + 1; Buffer[n].isSpace(); n++) ;
					if ((m = Buffer.indexOf(QString("°C").right(2), n+1)) != -1)
						TUnit->setCurrentIndex(0);
					else
					{
						TUnit->setCurrentIndex(1);
						if ((m = Buffer.indexOf("K", n+1)) == -1) m = Buffer.length() - 1;
					}
					while (Buffer[m].isSpace()) m--;
					Temp->setText(Buffer.mid(n, m-n+1));
					continue;
				}
				if (Buffer.indexOf("Laser frequency", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":", 0) + 1; Buffer[n].isSpace(); n++) ;
					LaserF->setText(Buffer.right(Buffer.length() - n));
					continue;
				}
				if ((n = Buffer.indexOf("Integration window from ", 0, Qt::CaseInsensitive)) != -1)
				{
					m = Buffer.indexOf("to", 0, Qt::CaseInsensitive);
					n+=24;
					WindowS->setText(Buffer.mid(n, m-n-1));
					WindowE->setText(Buffer.right(Buffer.length() - m - 3));
					continue;
				}
				if (Buffer.indexOf("Number of lines", 0, Qt::CaseInsensitive) != -1)
				{
					NumLines = Buffer.right(Buffer.length() - Buffer.indexOf(":", 0) - 1).toInt();
					if (Lines != 0) delete[] Lines;
					if (NumLines > 0) Lines = new SimLine[NumLines];
					else Lines = 0;
					continue;
				}
				break;
			case 1:
				if (l >= NumLines) continue;
				SList = Buffer.split("\t");
				if (SList.count() < 9) continue;
				for (n=0, Lines[l].Iso = 0; n < NI; n++) if (SList[0] == IsoN[n]) Lines[l].Iso = n; 
				Lines[l].Su = Lines[l].Sl = 0;
				for (n=0; n < Mol->getNumStates(); n++)
				{
					if ((St = Mol->getStateP(n))->getName() == SList[1]) Lines[l].Su = St;
					else if (St->getName() == SList[2]) Lines[l].Sl = St;
				}
				if (Lines[l].Su == 0 || Lines[l].Sl == 0) continue;
				Lines[l].vs = SList[3].toInt();
				Lines[l].Js = SList[4].toInt();
				Lines[l].vss = SList[5].toInt();
				Lines[l].Jss = SList[6].toInt();
				Lines[l].R = SList[7].toDouble();
				Lines[l++].I = SList[8].toDouble();
				break;
			case 2:
				SList = Buffer.split("\t");
				if (SList.count() >= 2) 
					Daten->AddValue(SList[0].toDouble(), SList[1].toDouble(), false);
				break;
		}
	}
	if (l < NumLines) l = NumLines;
	return true;
}

void SpectSimulation::RefreshMolecules()
{
	int n, N = MW->getNumMolecules();
	Molecule *MB;
	MolBox->clear();
	for (n=0; n<N; n++)
	{
		MolBox->addItem((MB = MW->getMolecule(n))->getName());
		if (MB == Mol) MolBox->setCurrentIndex(n);
	}
}

void SpectSimulation::SimBoxChanged(int nSim)
{
	if (Sim == nSim) return;
	if (Sim == 1 && nSim != 1) LaserF->setEnabled(false);
	else if (Sim != 1 && nSim == 1) LaserF->setEnabled(true);
	if (Sim == 2 && nSim != 2)
	{
		WindowE->setEnabled(false);
		WindowS->setEnabled(false);
	}
	if (Sim != 2 && nSim == 2)
	{
		WindowE->setEnabled(true);
		WindowS->setEnabled(true);
	}
	Sim = nSim;
}

void SpectSimulation::Simulate()
{
	printf("SpektSimulation::Simulate\n");
	//PaintScale(Bild->getPixmap(), Bild->contentsRect());
	double Rs = xStart->text().toDouble(), Re = xStop->text().toDouble(), RRes = 1.0 / XSF;
	double T = Temp->text().toDouble();
	printf("RRes=%e, XSF=%e\n", RRes, XSF);
	if (TUnit->currentIndex() == 0) T += C_K_C;
	double MRRes = 0.8 * Rs / C_c * sqrt(C_kB * T / (Mol->getIsoMass(0) * C_u));
	printf("T=%f, MRRes=%e\n", T, MRRes);
	if (RRes > MRRes) RRes = MRRes;
	int nP = int((Re - Rs) * XSF) + 1, i;
	double *Data = new double[nP], R, IRes = 0.1 / YSF;
	for (i=0; i < nP; i++) Data[i] = 0.0;
	if (Lines != 0) delete[] Lines;
	switch (Sim)
	{
		case 0:
			simulateAbsorption(Mol, Rs, Re, RRes, IRes, T, Data, NumLines, Lines);
			break;
		case 1:
			simulateLIF(Mol, Rs, Re, RRes, IRes, T, LaserF->text().toDouble(), Data, NumLines, Lines);
			break;
		case 2:
			simulateFilteredAbsorption(Mol, Rs, Re, RRes, IRes, T, WindowE->text().toDouble(), 
									   WindowS->text().toDouble(), Data, NumLines, Lines);
			break;
	}
	Daten->reinit();
	printf("Rs=%f, Re=%f, RRes=%f\n", Rs, Re, RRes);
	for (i=0, R=Rs; R <= Re; i++, R += RRes) Daten->AddValue(R, Data[i], false);
	delete[] Data;
	Paint();
	Changed();
	printf("Ende Simulate\n");
}

void SpectSimulation::simulateAbsorption(Molecule *Mol, double Rs, double Re, double RRes, double IRes,
										 double T, double *Data, int &NLines, SimLine *&Lines) const
{
	//printf("Simulate Absorption\n");
	bool QA = false;
    int vs, vss, Js, Jss, S, sJStep, JStep, JStartu, JStartl, n, m, I, NS, FMvs, FMvss, C, NumWFPoints = NumPoints;
	double R, f, DF = C_c * C_c / (-2.0 * C_kB * T), DFI, NA, RS = Rs - 1.0, RE = Re + 1.0, LI, HLF;
	ElState **ES;
	getStates(Mol, NS, ES);
	double ****TData[NS], ***FCF = 0, **ThPop, LRs, LRe, L_qsq = def_L_gsq, *LData;
	int NC[NS], MJ[NS], Mv[NS], NvP, NJP;
	TermTable *TT[NS];
	Potential *Pot[NS];
	FlexVector<SimLine> LineData(0);
	for (S=0; S < NS; S++)
	{
		Pot[S] = ES[S]->getPotential();
		if ((TT[S] = ES[S]->getTermTable()) == 0) 
			Pot[S]->calcTermEnergies(TT[S], cMaxv, cMaxJ, false);
		NC[S] = TT[S]->getNumComp();
		MJ[S] = TT[S]->getJMax();
		Mv[S] = TT[S]->getvMax();
		TData[S] = TT[S]->getData();
	}
	for (I = NLines = 0, IRes *= RRes; I < Iso->numIso; I++) 
	{
		NA = Iso->relNA[I];
		JStep = Mol->getJStep(I);
		JStartl = ES[0]->getJStart(I, 0);
		DFI = Mol->getIsoMass(I) * DF;
		TT[0]->getThermPopulation(ThPop, NvP, NJP, T, I, 0); 
		for (S=1; S < NS; S++) 
		{
			JStartu = ES[S]->getJStart(I, 0);
			if (NC[S] >= 2 ? (Js = ES[S]->getJStart(I, 1)) < JStartu : false) JStartu = Js;
			if (ES[S]->getLambda() > 0 && NC[S] > 1) 
			{
				QA = true;
				sJStep = 1;
			}
			else 
			{
				QA = false;
				sJStep = JStep;
			}
			for (Js = JStartu; Js <= MJ[S]; Js += sJStep) for (Jss = Js - 1; Jss <= Js + 1; Jss++)
					if (Jss >= JStartl && (QA || Jss != Js) && Jss <= MJ[0])
			{
				//printf("I=%d, S=%d, Js=%d\n", I, S, Js);
				if (FCF != 0) 
				{
					Destroy(FCF, FMvs, FMvss);
					FCF = 0;
				}
				if (Js == Jss) 
				{
					C=1;
					HLF = 2.0;
				}
				else 
				{
					C = 0;
					HLF = 1.0;
				}
				FMvs = Mv[S];
				FMvss = Mv[0];
				for (vs=n=0, R = Rs; (vs <= FMvs ? TData[S][C][I][vs][Js] != 0 : false); vs++) 
					for (vss=0; (vss <= FMvss ? TData[0][0][I][vss][Jss] != 0 : false); vss++) 
				{
					f = (FCF == 0 ? TData[S][C][I][vs][Js] - TData[0][0][I][vss][Jss] 
						          : FCF[vs][vss][1]);
					if (f < RS || f > RE) continue;
                    if (FCF == 0) Pot[0]->getFCF(Pot[S], I, Js, Jss, FMvs, FMvss, FCF, NumWFPoints);
					if ((LI = NA * HLF * ThPop[vss][Jss] * FCF[vs][vss][0]) < IRes) continue;
					while (R<f) 
					{
						R += RRes;
						n++;
					}
					while (R>f) 
					{
						R -= RRes;
						n--;
					}
					voigt(f, LI, LRs = R, LRe, LData, RRes, IRes, DFI / (f*f), L_qsq);
					while (R > LRs)
					{
						R -= RRes;
						n--;
					}
					if (n<0)
					{
						R = Rs;
						m = -n;
						n=0;
					}
					else m=0;
					while (R <= LRe && R <= Re)
					{
						Data[n++] += LData[m++];
						R += RRes;
					}
					delete[] LData;
					LineData[NLines].I = LI;
					LineData[NLines].Iso = I;
					LineData[NLines].Js = Js;
					LineData[NLines].Jss = Jss;
					LineData[NLines].R = f;
					LineData[NLines].Sl = ES[0];
					LineData[NLines].Su = ES[S];
					LineData[NLines].vs = vs;
					LineData[NLines++].vss = vss;
				}
			}
		}
		Destroy(ThPop, NvP);
	}
	if (FCF != 0) Destroy(FCF, FMvs, FMvss);
	delete[] ES;
	Lines = LineData.getData();
}

void SpectSimulation::simulateFilteredAbsorption(Molecule *Mol, double Rs, double Re, double RRes,
		 										 double IRes, double T, double Ws, double We, 
			  									 double *Data, int &NLines, SimLine *&Lines) const
{
	printf("SimulateFilteredAbsorption\n");
	bool QA = false;
    int v, J, vs, vss, Js, Jss, S, sJStep, JStep, JStartu, JStartl, n, m, I, NS, FMvs, FMvss[3], C, cn, NumWFPoints = NumPoints;
	double R, f, DF = C_c * C_c / (-2 * C_kB * T), NA, RS = Rs - 1.0, RE = Re + 1.0, LI, HLF, M = 0.0;
	ElState **ES; 
	getStates(Mol, NS, ES);
	double ****TData[NS], ***FCF[3], **ThPop, LRs, LRe, L_gsq = def_L_gsq, *LData, *EProf, PRs, PRe;
	double WS = Ws - 1.0, WE = We + 1.0, E, LF = 0.0, IS, CR;
	int NC[NS], MJ[NS], Mv[NS], NvP, NJP, Fie, Fif, i, d;
	TermTable *TT[NS];
	Potential *Pot[NS];
	FlexVector<SimLine> LineData(0);
	for (n=0; n<3; n++) FCF[n] = 0;
	for (S=0; S < NS; S++)
	{
		Pot[S] = ES[S]->getPotential();
		if ((TT[S] = ES[S]->getTermTable()) == 0) 
			Pot[S]->calcTermEnergies(TT[S], cMaxv, cMaxJ, false);
		NC[S] = TT[S]->getNumComp();
		MJ[S] = TT[S]->getJMax();
		Mv[S] = TT[S]->getvMax();
		TData[S] = TT[S]->getData();
	}
	int FLines[2*Mv[0]][3], NFL;
	for (I = NLines = 0, IRes *= RRes; I < Iso->numIso; I++) 
	{
		NA = Iso->relNA[I];
		JStep = Mol->getJStep(I);
		JStartl = ES[0]->getJStart(I, 0);
		//DFI = Mol->getIsoMass(I) * DF;
		TT[0]->getThermPopulation(ThPop, NvP, NJP, T, I, 0); 
		for (S=1; S < NS; S++) 
		{
			JStartu = ES[S]->getJStart(I, 0);
			if (NC[S] >= 2 ? (Js = ES[S]->getJStart(I, 1)) < JStartu : false) JStartu = Js;
			if (ES[S]->getLambda() > 0 && NC[S] > 1) 
			{
				QA = true;
				sJStep = 1;
			}
			else 
			{
				QA = false;
				sJStep = JStep;
			}
			for (Js = JStartu; Js <= MJ[S]; Js += sJStep) 
			{
				for (n=0; n<3; n++) 
				{
					if (FCF[n] != 0) 
					{
						Destroy(FCF[n], FMvs, FMvss[n]);
						FCF[n] = 0;
					}
					FMvss[n] = Mv[0];
				}
				FMvs = Mv[S];
				for (Jss = Js - 1, Fie = 0; Jss <= Js + 1; Jss++, Fie++)
					if (Jss >= JStartl && (QA || Jss != Js) && Jss <= MJ[0])
				{
					if (Js == Jss) 
					{
						C=1;
						HLF = 2.0;
					}
					else 
					{
						C = 0;
						HLF = 1.0;
					}
					for (vs=n=0, R = Rs; (vs <= FMvs ? TData[S][C][I][vs][Js] != 0 : false); vs++) 
						for (vss=0; (vss <= FMvss[Fie] ? TData[0][0][I][vss][Jss] != 0 : false); vss++)
					{
						f = (FCF[Fie] == 0 ? TData[S][C][I][vs][Js] - TData[0][0][I][vss][Jss] 
							          : FCF[Fie][vs][vss][1]);
						if (f < RS || f > RE 
									|| (LI = NA * HLF * ThPop[vss][Jss] * FCF[Fie][vs][vss][0]) < IRes)
							continue;
						for (NFL = 0, J = Js -1 + (Fif = (Js == Jss ? 1 : 0)); 
											   Fif <= 2; J+=2, Fif += 2) 
							for (v=0; v <= FMvss[Fif]; v++) 
						{
							E = (FCF[Fif] == 0 ? TData[S][C][I][vs][Js] - TData[0][0][I][v][J]
								               : FCF[Fif][vs][v][1]);
							if (E > WS && E < WE) 
							{
								FLines[NFL][0] = v;
								FLines[NFL][1] = J;
								FLines[NFL++][2] = Fif;
							}
						}
						if (NFL == 0) continue;
						if (FCF[Fie] == 0) 
                            Pot[0]->getFCF(Pot[S], I, Js, Jss, FMvs, FMvss[Fie], FCF[Fie], NumWFPoints);
						for (m=0, E=0.0; m < NFL; m++) 
						{
							if (FCF[FLines[m][2]] == 0)
								Pot[0]->getFCF(Pot[S], I, Js, FLines[m][1], FMvs, FMvss[FLines[m][2]],
                                               FCF[FLines[m][2]], NumWFPoints);
							E += FCF[FLines[m][2]][vs][FLines[m][0]][1];
						}
						if (E * LI < IRes) continue;
						while (LF < f) 
						{
							LF += RRes;
							n++;
						}
						while (LF > f) 
						{
							LF -= RRes;
							n--;
						}
						for (cn = n, CR = LF, d=1, IS = IRes; true; n+=d, LF += (d==1 ? RRes : -RRes))
						{
							if (IS < IRes) 
							{
								if (d==1) 
								{
									d=-1;
									n = cn - 1;
									LF = CR - RRes;
								}
								else break;
							}
							cExcProfile(f, LI, LF, PRs, PRe, EProf, RRes, IRes, DF / (f*f), L_gsq);
							for (i=0, IS = 0.0; i < NFL; i++)
							{
								flLine(FCF[FLines[i][2]][vs][FLines[i][0]][1], f,
									   FCF[FLines[i][2]][vs][FLines[i][0]][0], 
									   LRs = FCF[FLines[i][2]][vs][FLines[i][0]][1] + f - PRs, LRe,
			 						   LData, RRes, IRes, L_gsq, PRs, PRe, EProf);
								for (R = LRs, m=0; R < Ws; R += RRes, m++) ;
								if (R > We)
								{
									m--;
									R -= RRes;
								}
								while (R <= LRe && R <= We)
								{
									IS += LData[m++];
									R += RRes;
								}
								delete[] LData;
							}
							if (n == cn) M = IS;
							else if (n == cn - 1 && M < IS) M = IS;
							delete[] EProf;
							Data[n] += IS;
						}
						LineData[NLines].I = M;
						LineData[NLines].Iso = I;
						LineData[NLines].Js = Js;
						LineData[NLines].Jss = Jss;
						LineData[NLines].R = f;
						LineData[NLines].Sl = ES[0];
						LineData[NLines].Su = ES[S];
						LineData[NLines].vs = vs;
						LineData[NLines++].vss = vss;
					}
				}
			}
		}
		Destroy(ThPop, NvP);
	}
	for (n=0; n<3; n++) if (FCF[n] != 0) Destroy(FCF[n], FMvs, FMvss[n]);
	delete[] ES;
	Lines = LineData.getData();
}

void SpectSimulation::simulateLIF(Molecule *Mol, double Rs, double Re, double RRes, double IRes, double T, 
					              double LF, double *Data, int &NLines, SimLine *&Lines) const
{
	//printf("SimulateLIF\n");
	bool QA = false;
    int v, J, vs, vss, Js, Jss, S, sJStep, JStep, JStartu, JStartl, n, m, I, NS, FNvs, FNvss[3], C, NumWFPoints = NumPoints;
	double R, f, DF = C_c * C_c / (-2 * C_kB * T), NA, RS = LF - 1.0, RE = LF + 1.0, LI, HLF, M;
	ElState **ES ;
	getStates(Mol, NS, ES);
	double ****TData[NS], ***FCF[3], **ThPop, LRs, LRe, L_gsq = def_L_gsq, *LData, *EProf, PRs, PRe;
	int NC[NS], MJ[NS], Mv[NS], NvP, NJP, Fie, Fif, SRs, SRe;
	TermTable *TT[NS];
	Potential *Pot[NS];
	FlexVector<SimLine> LineData(0);
	for (n=0; n<3; n++) FCF[n] = 0;
	for (S=0; S < NS; S++)
	{
		Pot[S] = ES[S]->getPotential();
		if ((TT[S] = ES[S]->getTermTable()) == 0) 
			Pot[S]->calcTermEnergies(TT[S], cMaxv, cMaxJ, false);
		NC[S] = TT[S]->getNumComp();
		MJ[S] = TT[S]->getJMax();
		Mv[S] = TT[S]->getvMax();
		TData[S] = TT[S]->getData();
	}
	for (I = NLines = 0, IRes *= RRes; I < Iso->numIso; I++) 
	{
		NA = Iso->relNA[I];
		JStep = Mol->getJStep(I);
		JStartl = ES[0]->getJStart(I, 0);
		//DFI = Mol->getIsoMass(I) * DF;
		TT[0]->getThermPopulation(ThPop, NvP, NJP, T, I, 0); 
		for (S=1; S < NS; S++) 
		{
			JStartu = ES[S]->getJStart(I, 0);
			if (NC[S] >= 2 ? (Js = ES[S]->getJStart(I, 1)) < JStartu : false) JStartu = Js;
			if (ES[S]->getLambda() > 0 && NC[S] > 1) 
			{
				QA = true;
				sJStep = 1;
			}
			else 
			{
				QA = false;
				sJStep = JStep;
			}
			for (Js = JStartu; Js <= MJ[S]; Js += sJStep) 
			{
				for (n=0; n<3; n++) if (FCF[n] != 0) 
				{
					Destroy(FCF[n], FNvs, FNvss[n]);
					FCF[n] = 0;
				}
				for (Jss = Js - 1, Fie = 0; Jss <= Js + 1; Jss++, Fie++)
					if (Jss >= JStartl && (QA || Jss != Js) && Jss <= MJ[0])
				{
					if (Js == Jss) 
					{
						C=1;
						HLF = 2.0;
					}
					else 
					{
						C = 0;
						HLF = 1.0;
					}
					if (FCF[Fie] == 0)
					{
						FNvs = Mv[S] + 1;
						FNvss[Fie] = Mv[0] + 1;
					}
					for (vs=n=0, R = Rs; (vs < FNvs ? TData[S][C][I][vs][Js] != 0 : false); vs++) 
						for (vss=0; (vss < FNvss[Fie] ? TData[0][0][I][vss][Jss] != 0 : false); vss++)
					{
						f = (FCF[Fie] == 0 ? TData[S][C][I][vs][Js] - TData[0][0][I][vss][Jss] 
							          : FCF[Fie][vs][vss][1]);
						if (f < RS || f > RE) continue;
						if (FCF[Fie] == 0) 
                            Pot[0]->getFCF(Pot[S], I, Js, Jss, FNvs, FNvss[Fie], FCF[Fie], NumWFPoints);
						if ((LI = NA * HLF * ThPop[vss][Jss] * FCF[Fie][vs][vss][0]) < IRes) continue;
						cExcProfile(f, LI, LF, PRs, PRe, EProf, RRes, IRes, DF / (f*f), L_gsq);
						for (R = PRs, m=0, M=0.0; R <= PRe && (R <= f || R <= LF); m++, R += RRes) 
							if (EProf[m] > M) M = EProf[m];
						if (M < IRes) continue;
						SRs = Rs + f - PRs;
						SRe = Re - f + PRe;
						if (Js != Jss)
						{
							if (Fie == 0) 
							{
								Fif = 2;
								J = Js + 1;
							}
							else 
							{
								Fif = 0;
								J = Js - 1;
							}
							if (FCF[Fif] == 0 && J >= JStartl && J <= MJ[S]) 
                                Pot[0]->getFCF(Pot[S], I, Js, J, FNvs, FNvss[Fif], FCF[Fif], NumWFPoints);
						}
						for (J = Js - 1 + (Fif = (Js == Jss ? 1 : 0)); Fif <= 2; J+=2, Fif += 2) 
							if (J >= JStartl && J <= MJ[S]) for (v=0; v < FNvss[Fif]; v++) 
									if (FCF[Fif][vs][v][1] > SRs && FCF[Fif][vs][v][1] < SRe 
																   && M * FCF[Fif][vs][v][0] > IRes)
						{
							LRs = FCF[Fif][vs][v][1] + f - PRs;
							while (R < LRs) 
							{
								R += RRes;
								n++;
							}
							while (R > LRs) 
							{
								R -= RRes;
								n--;
							}
							flLine(FCF[Fif][vs][v][1], f, FCF[Fif][vs][v][0], LRs = R, LRe, LData,
								   RRes, IRes, L_gsq, PRs, PRe, EProf);
							while (R > LRs)
							{
								R -= RRes;
								n--;
							}
							if (n<0)
							{
								R = Rs;
								m = -n;
								n=0;
							}
							else m=0;
							while (R <= LRe && R <= Re)
							{
								Data[n++] += LData[m++];
								R += RRes;
							}
							delete[] LData;
							LineData[NLines].I = M * FCF[Fif][vs][v][0];
							LineData[NLines].Iso = I;
							LineData[NLines].Js = Js;
							LineData[NLines].Jss = J;
							LineData[NLines].R = FCF[Fif][vs][v][1];
							LineData[NLines].Sl = ES[0];
							LineData[NLines].Su = ES[S];
							LineData[NLines].vs = vs;
							LineData[NLines++].vss = v;
						}
						delete[] EProf;
					}
				}
			}
		}
		Destroy(ThPop, NvP);
	}
	for (n=0; n<3; n++) if (FCF[n] != 0) Destroy(FCF[n], FNvs, FNvss[n]);
	delete[] ES;
	Lines = LineData.getData();
}

void SpectSimulation::voigt(double f, double I, double &Rs, double &Re, double *&Data, double RRes,
							double IRes, double DF, double L_gsq) const
{
	double D, V, wD, wL, *DP, *LP;
	int i, j, k, nD, nL;
	DoubleVector R;
	for (i=0, wD = Rs - f, D = IRes; D >= IRes; i++, wD += RRes)
		R[i] = (D = I * exp(DF * wD * wD));
	Re = f + wD - RRes;
	for (j=-1, wD = Rs - f, D = IRes; D >= IRes; j--, wD -= RRes)
		R[j] = (D = I * exp(DF * wD * wD));
	Rs = f - wD + RRes;
	DP = R.getData();
	nD = i-j-1;
	R.clear();
	R[0] = 1.0 / sqrt(L_gsq);
	for (j=1, k=-1, wL = RRes, V = D = IRes / I; V >= D; j++, k--, wL += RRes) 
	{
		R[j] = V = 1.0 / sqrt(wL * wL + L_gsq);
		R[k] = V;
	}
	Rs -= (wL -= RRes);
	Re += wL;
	LP = R.getData();
	nL = 2*j-1;
	R.clear();
	for (i=0; i < nD; i++) for (j=0; j < nL; j++) R[i+j] += DP[i] * LP[j];
	Data = R.getData();
	delete[] DP;
	delete[] LP;
}

double SpectSimulation::voigtP(double P, double f, double I, double RRes, double IRes, double DF,
							   double L_gsq) const
{
	double D, V, wD, wL, *DP, *LP, Rs;
	int i, j, k=0, nD, nL;
	DoubleVector R;
	for (Rs = P; Rs < f; Rs += RRes, k++) ;
	while (Rs > f) 
	{
		Rs -= RRes;
		k--;
	}
	for (i=0, wD = Rs - f, D = IRes; D >= IRes; i++, wD += RRes)
		R[i] = (D = I * exp(DF * wD * wD));
	for (j=-1, wD = Rs - f, D = IRes; D >= IRes; j--, wD -= RRes)
		R[j] = (D = I * exp(DF * wD * wD));
	DP = R.getData();
	k = (nD = i-j-1) - k;
	R.clear();
	R[0] = 1.0 / sqrt(L_gsq);
	for (j=1, k=-1, wL = RRes, V = D = IRes / I; V >= D; j++, k--, wL += RRes) 
	{
		R[j] = V = 1.0 / sqrt(wL * wL + L_gsq);
		R[k] = V;
	}
	LP = R.getData();
	nL = 2*j-1;
	if ((i=k-j-1) < 0)
	{
		j = -i;
		i = 0;
	}
	else j=0;
	for (V=0.0; i < nD && j < nL; i++, j++) V += DP[i] * LP[j];
	delete[] DP;
	delete[] LP;
	return V;
}

bool SpectSimulation::writeData(QString FileName)
{
	QFile Datei(FileName);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
	int i;
	double *P;
	S << "Simulated spectrum:\n";
	S << "Type: " << SimBox->currentText() << "\n";
	S << "Energy range from " << xStart->text() << " to " << xStop->text() << "\n";
	S << "Intensity range form " << yStart->text() << " to " << yStop->text() << "\n";
	S << "Molecule: " << MolBox->currentText();
	S << "Temperature: " << Temp->text() << " " << TUnit->currentText() << "\n";
	if (SimBox->currentIndex() == 1) S << "Laser frequency: " << LaserF->text() << "\n";
	else if (SimBox->currentIndex() == 2) 
		S << "Integration window from " << WindowS->text() << " to " << WindowE->text() << "\n";
	S << "Number of lines: " << NumLines << "\n";
	S << "Number of points: " << Daten->GetDSL() << "\n\n";
	S << "Table of lines:\n";
	S << "Molecule:\tupper state:\tlower state:\tJ'\tv'\tJ''\tv''\tfrequency:\tintensity:\n";
	for (i=0; i < NumLines; i++) 
			S << Iso->getIsoName(Lines[i].I) << "\t" << Lines[i].Su->getName() << "\t"
				<< Lines[i].Sl->getName() << "\t" << Lines[i].vs << "\t" << Lines[i].Js << "\t"
				<< Lines[i].vss << "\t" << Lines[i].Jss << "\t" << QString::number(Lines[i].R, 'f', 11)
				<< "\t" << QString::number(Lines[i].I, 'g', 8) << "\n";
	S << "\n" << "Table of data points:\n";
	S << "R [A]:\t<<intensity [arbitrary units]:\n";
	for (i=0; i < Daten->GetDSL(); i++)
	{
		P = Daten->GetPoint(i);
		S << QString::number(P[0], 'f', 11) << "\t" << QString::number(P[1], 'g', 8) << "\n";
	}
	return true;
}
