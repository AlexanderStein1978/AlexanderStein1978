//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "addspektrum.h"
#include "markedpeak.h"
#include "ResidualFit.h"
#include "termenergy.h"
#include "fitdata.h"
#include "adddialog.h"
#include "elstate.h"
#include "utils.h"
#include "termtable.h"
#include "perturbation.h"
#include "point.h"
#include "tableline.h"
#include "molecule.h"
#include "isotab.h"
#include "doubleaddspektrum.h"
#include "fit.h"

#include <cmath>

#include <QMessageBox>


AddSpectrum::AddSpectrum(MainWindow* MW, AddDialog *D) : DiagWindow(AddSpect, MW), m_contrast(0), Data(0), Ri(0.0), Ra(0.0), Res(0.0), PTData(0), OTData(0), UEn(0),
    m_IntensF(0), m_ownResidualFit(false)
{
	ResFit = 0;
	DAddSpect = 0;
	Dialog = D;
	FitD = 0;
	PFitD = 0;
	marked = new MarkedPeak[1000];
	NMarked = type = TDC = Comp = v = Iso = NP = NJ = JStart = JStep = NC = 0;
	NPert = NOC = NOIso = NOJ = NOv = 0;
	NCP = NvP = NJP = EShift = 0;
	TT = OffTT = 0;
	UCP = 0;
	State = 0;
	Source = 0;
	vAss = 0;
	CAss = 0;
	DAddSpect = 0;
	NMBAUL = -1;
    CompBox->addItem("No shift");
	PMenu = new QMenu();
	markLevelsAct = new QAction("&Mark levels", this);
	removeLevelsAct = new QAction("&Remove levels", this);
	updateLevelsAct = new QAction("&Update levels", this);
	hideLevelsAct = new QAction("&Hide levels", this);
	showLevelsAct = new QAction("&Show levels", this);
	autoUpdateAct = new QAction("&Auto update", this);
	markLevelsAct->setStatusTip("Mark the rows in the fit dataset which belong to the selected levels"); 
	removeLevelsAct->setStatusTip("Delete the selected levels from the fit dataset");
	updateLevelsAct->setStatusTip(
		"Update the energy in the fit dataset of all / the selected levels");
	hideLevelsAct->setStatusTip("Hide the selected levels or all levels if no levels are selected");
	showLevelsAct->setStatusTip("Show all hidden levels");
	autoUpdateAct->setStatusTip(
		"Generate a suggestion for the update of the assigments of the whole\
			vibrational state");
	connect(markLevelsAct, SIGNAL(triggered()), this, SLOT(markLevels()));
	connect(removeLevelsAct, SIGNAL(triggered()), this, SLOT(removeLevels()));
	connect(updateLevelsAct, SIGNAL(triggered()), this, SLOT(updateLevels()));
	connect(hideLevelsAct, SIGNAL(triggered()), this, SLOT(hideLevels()));
	connect(showLevelsAct, SIGNAL(triggered()), this, SLOT(showLevels()));
	connect(autoUpdateAct, SIGNAL(triggered()), this, SLOT(autoUpdateLevels()));
	connect(this, SIGNAL(clicked(double,double)), this, SLOT(ImageClicked(double,double)));
	PMenu->addAction(markLevelsAct);
	PMenu->addAction(hideLevelsAct);
	PMenu->addAction(showLevelsAct);
	PMenu->addAction(updateLevelsAct);
	PMenu->addAction(removeLevelsAct);
	PMenu->addSeparator();
	PMenu->addAction(autoUpdateAct);
	AutoUpdateMenu = new QMenu();
	QAction *acceptResultAct = new QAction("&Accept result", this);
	QAction *rejectResultAct = new QAction("&Reject result", this);
	acceptResultAct->setStatusTip(
		"Apply the suggested assignment change and save the changed level assignment to disk");
	rejectResultAct->setStatusTip(
		"Withdraw the suggestion and do not change anything to the level assignment");
	connect(acceptResultAct, SIGNAL(triggered()), this, SLOT(acceptAutoUpdateResult()));
	connect(rejectResultAct, SIGNAL(triggered()), this, SLOT(rejectAutoUpdateResult()));
	AutoUpdateMenu->addAction(acceptResultAct);
	AutoUpdateMenu->addAction(rejectResultAct);
	connect(Bild, SIGNAL(RightClicked(QPoint*)), this, SLOT(showPMenu(QPoint*)));
	connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(mark(QRect*)));
	connect(SourceBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setEnergyReference()));
	connect(CompBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setEnergyReference()));
    connect(m_contrastSlider, SIGNAL(valueChanged(int)), this, SLOT(contrastChanged(int)));
    connect(m_intensitySlider, SIGNAL(valueChanged(int)), this, SLOT(Paint()));
}

AddSpectrum::~AddSpectrum()
{
	delete marked;
	if (EShift != 0) delete[] EShift;
	if (PFitD != 0)
	{
		delete[] PFitD;
		delete[] OTData;
		delete[] NCP;
		delete[] NvP;
		delete[] NJP;
	}
	delete PMenu;
	if (UCP != 0)
	{
		int n;
		for (n=0; n < NPert; n++) delete[] UCP[n];
		delete[] UCP;
	}
	if (ResFit != 0 && m_ownResidualFit) delete ResFit;
    if (m_IntensF != 0) delete[] m_IntensF;
}

void AddSpectrum::AssignmentsAccepted(FitData *CFD)
{
	int n, f=-2;
	if (CFD == FitD) f=-1;
	else for (n=0; n < NPert; n++) if (CFD == PFitD[n]) f=n;
	if (f==-2) return;
	for (n=0; n < NMarked; n++) 
		if ((marked[n].LevelState == 0 || marked[n].LevelState == 2) 
			&& marked[n].StateNum == f) marked[n].LevelState = 1;
	Paint();
}

void AddSpectrum::calcEJ(double& rE, int& rJ, int J, int P)
{
	double x, y;
	rJ = JStart + JStep * J;
	x = UEn[J] - Ri + double(P) * Res;
	if (P>0 && P < NP - 1) 
		ParabInterpol(x - Res, Data[J][P-1], x, Data[J][P], x + Res, Data[J][P+1], rE, y);
	else rE = x;
}

void AddSpectrum::contrastChanged(int i_c)
{
    double stretch = pow(2.0, i_c - m_contrast);
	int cf = static_cast<int>(pow(2.0, i_c) + 0.1);
    m_intensitySlider->blockSignals(true);
    m_intensitySlider->setMaximum(256);
    m_intensitySlider->setMinimum(-cf * 256);
    m_intensitySlider->setValue(static_cast<int>(stretch * static_cast<double>(m_intensitySlider->value() + 128)) - 128);
    m_intensitySlider->blockSignals(false);
    m_contrast = i_c;
    Paint();
}

void AddSpectrum::getMarked(int &N, MarkedPeak *&MarkedLevels)
{
	N = NMarked;
	MarkedLevels = marked;
}

void AddSpectrum::hideLevels()
{
	if (NMarked == 0) return;
	int n;
	for (n=0; (n < NMarked ? marked[n].LevelState < 2 : false); n++) ;
	if (n == NMarked) 
	{
		for (n=0; n < NMarked; n++) if (marked[n].LevelState == 1) 
			marked[n].LevelState = -1;
	}
	else for (n=0; n < NMarked; n++)
	{
		if (marked[n].LevelState == 2) marked[n].LevelState = 0;
		else if (marked[n].LevelState == 3) marked[n].LevelState = -1;
	}
	Paint();
}

void AddSpectrum::fillTEStruct(TermEnergy& TE, int &mt, int J, int M)
{
	int n, mE, c, i, mc = -1, mv = -1;
	calcEJ(TE.E, TE.J, J, M);
	if (OffTT != 0) 
	{
		for (n=0, mt = -1, mE = fabs(TE.E - OTData[Comp][Iso][v][TE.J]); n < NPert; n++)
			if (TE.J < NJP[n]) for (c=0; c < NCP[n]; c++) if (UCP[n][c]) 
				for (i=0; i < NvP[n]; i++)
					if (fabs(TE.E - PTData[n][c][Iso][i][TE.J]) < mE 
						&& PTData[n][c][Iso][i][TE.J] != 0.0)
		{
			mE = fabs(TE.E - PTData[n][c][Iso][i][TE.J]);
			mc = c;
			mv = i;
			mt = n;
		}
	}
	else mt = -1;
	if (mt == -1 ? FitD != 0 : PFitD[mt] != 0)
	{
		TE.Iso = Iso;
		TE.v = (mt == -1 ? (vAss != 0 ? vAss[J] : v) : mv);
		TE.err = 0.0201;
		TE.dev = 0.0;
		TE.DevR = 0.0;
		TE.ef = ((mt == -1 ? Comp == 0 : mc==0) ? true : false);
		TE.FC = (CAss != 0 ? CAss[J] : -1);
        TE.State = State;
	}
}

void AddSpectrum::addLevel(TermEnergy& TE, int& mt, int J, int M)
{
	marked[NMarked].JPos = J;
	marked[NMarked].EPos = M;
	marked[NMarked].LevelState = 2;
	if (mt == -1) 
	{
		if (FitD != 0) marked[NMarked].Row = FitD->addMarkedLevel(TE, Source);
	}
	else if (PFitD[mt] != 0) marked[NMarked].Row = PFitD[mt]->addMarkedLevel(TE, Source);
	marked[NMarked].CalcFromProgression = false;
	marked[NMarked++].StateNum = mt;
	if (NMarked == 1) Dialog->enableOptimize(true);
}

void AddSpectrum::ImageClicked(double dJ, double E)
{
	if (NJ == 0) return;
	int J = (JStep == 2 ? int(0.5 * (dJ - double(JStart - 1))) : int(dJ) - JStart), n, M;
	if (J < 0 || J >= NJ) return;
	int P = int((UEn != 0 ? E + Ri : E - Ri) / Res) - EShift[J], mt, i;
	if (type == Absorption)
	{
		if (P>0 ? Data[J][P-1] < Data[J][P] : Data[J][P] < Data[J][P+1])
			for (M=P; (M>0 ? Data[J][M-1] < Data[J][M] : false); M--) ;
		else for (M=P; (M < NP - 1 ? Data[J][M+1] < Data[J][M] : false) ; M++) ;
	}
	else
	{
		if (P>0 ? Data[J][P-1] <= Data[J][P] : Data[J][P] < Data[J][P+1])
			for (M=P; (M < NP  - 1 ? Data[J][M+1] > Data[J][M] : false); M++) ;
		else for (M=P; (M>0 ? Data[J][M-1] > Data[J][M] : false); M--) ;
	}
	TermEnergy TE;
	fillTEStruct(TE, mt, J, M);
	for (n=0; (n < NMarked ? marked[n].JPos != J || marked[n].EPos != M : false); n++) ;
	if (n < NMarked)
	{
		switch (marked[n].LevelState)
		{
			case -1:
				marked[n].LevelState = 3;
				break;
			case 0:
				marked[n].LevelState = 2;
				break;
			case 1:
				marked[n].LevelState = 3;
				break;
			case 2:
				if (marked[n].StateNum == -1 && FitD != 0) 
					FitD->removeMarkedLevel(TE, Source);
				else if (marked[n].StateNum >= 0 ? PFitD[marked[n].StateNum] != 0 : false)
					PFitD[marked[n].StateNum]->removeMarkedLevel(TE, Source);
				for (i=0; i < NMarked; i++) 
					if (marked[i].StateNum == marked[n].StateNum 
						&& marked[i].Row > marked[n].Row) marked[i].Row--;
				for (n++; n < NMarked; n++)
				{
					marked[n-1].JPos = marked[n].JPos;
					marked[n-1].EPos = marked[n].EPos;
					marked[n-1].LevelState = marked[n].LevelState;
					marked[n-1].Row = marked[n].Row;
					marked[n-1].StateNum = marked[n].StateNum;
					marked[n-1].CalcFromProgression = marked[n].CalcFromProgression;
				}
				NMarked--;
				if (NMarked == 0) Dialog->enableOptimize(false);
				break;
			case 3:
				marked[n].LevelState = 1;
				break;
		}
	}
	else if (NMarked < 1000) addLevel(TE, mt, J, M);
	Paint();
	if (mt == -1 ? FitD != 0 : PFitD[mt] != 0)
	{
		if (mt == -1)
		{
			if (!FitD->isVisible()) FitD->show();
			FitD->activateWindow();
			FitD->setFocus();
		}
		else
		{
			if (!PFitD[mt]->isVisible()) PFitD[mt]->show();
			PFitD[mt]->activateWindow();
			PFitD[mt]->setFocus();
		}
		activateWindow();
		setFocus();
	}
}

void AddSpectrum::KeyPressEvent(QKeyEvent* K)
{
    if (K->key() == Qt::Key_H)
	{
		hideLevels();
		K->accept();
	}
	else if (K->key() == Qt::Key_S)
	{
		showLevels();
		K->accept();
	}
	else DiagWindow::KeyPressEvent(K);
}

void AddSpectrum::mark(QRect* R)
{
	if (ZoomB->isChecked()) return;
	int Em, EM, Jm, JM, n;
	if (UEn != 0)
	{
		Em = int((YO - double(R->bottom()) + Ri * YSF) / (YSF * Res));
		EM = int((YO - double(R->top()) + Ri * YSF) / (YSF * Res));
	}
	else
	{
		Em = int((-Ri * YSF - double(R->bottom()) + YO) / (YSF * Res));
		EM = int((-Ri * YSF - double(R->top()) + YO) / (YSF * Res));
	}
	if (JStep == 2)
	{
		Jm = int(0.5 * ((double(R->left()) - XO) / XSF - double(JStart - 1)));
		JM = int(0.5 * ((double(R->right()) - XO) / XSF - double(JStart - 1)));
	}
	else
	{
		Jm = int((double(R->left()) - XO) / XSF) - JStart;
		JM = int((double(R->right()) - XO) / XSF) -JStart;
	}
	for (n=0; n < NMarked; n++)
	{
		if (marked[n].JPos >= Jm && marked[n].JPos <= JM 
			&& marked[n].EPos >= Em - EShift[marked[n].JPos] 
			&& marked[n].EPos <= EM - EShift[marked[n].JPos])
		{
			if (marked[n].LevelState == -1 || marked[n].LevelState == 1) 
				marked[n].LevelState = 3;
			else if (marked[n].LevelState == 0) marked[n].LevelState = 2;
		}
		else if (marked[n].LevelState == 3) marked[n].LevelState = 1;
	}
	Paint();
}

void AddSpectrum::markLevels()
{
	int n, N[NPert + 1], *RNum[NPert + 1];
	if (FitD == 0) FitD = State->getFitData();
	for (n=0; n <= NPert; n++) N[n] = 0;
	for (n=0; n < NMarked; n++) if (marked[n].LevelState >= 2) 
		N[marked[n].StateNum + 1]++;
	for (n=0; n <= NPert; n++) if (N[n] > 0)
	{
		RNum[n] = new int[N[n]];
		N[n] = 0;
	}
	for (n=0; n < NMarked; n++) if (marked[n].LevelState >= 2) 
		RNum[marked[n].StateNum + 1][N[marked[n].StateNum + 1]++] = marked[n].Row;
	if (FitD != 0 && N[0] > 0) FitD->MarkLines(RNum[0], N[0]);
	for (n=0; n < NPert; n++) if (PFitD[n] != 0 && N[n+1] > 0) 
		PFitD[n]->MarkLines(RNum[n+1], N[n+1]);
	for (n=0; n <= NPert; n++) if (N[n] > 0) delete[] RNum[n];
}

void AddSpectrum::PSpektrum(QPainter& P, const QRect& A, bool PrintFN)
{
	if (XMax == XMin || YMax == YMin) return;
	QRgb *PP;
    double /*ISc,*/ mI, MI, JISc, curVal, CF = pow(2.0, m_contrast);
	int i, b, m1, m2, n, m, J, **Pixel = CreateInt(NJ, NP), intensityOffset = m_intensitySlider->value();
	QImage *Pict = new QImage(NJ, NP, QImage::Format_RGB32);
	double xmin = xStart->text().toDouble(), xmax = xStop->text().toDouble();
	double ymin = yStart->text().toDouble(), ymax = yStop->text().toDouble();
	double xsc = NJ / (XMax - XMin);
	double ysc = NP / (YMax - YMin);
	int JB = int((xmin - XMin) * xsc), EB = int((ymin - YMin) * ysc), EBJ, ESJ;
    int JS = JB + int((xmax - xmin) * xsc) + 1, ES = EB + int((ymax - ymin) * ysc) + 1, EC = NP / 2, *BFJ = new int[JS+1];
	if (EB < 0) EB = 0;
	if (JB < 0) JB = 0;
	if (ES >= NP) ES = NP - 1;
	if (JS >= NJ) JS = NJ - 1;
    /*for (mI = 1e9, MI = -1e9, J = JStart + JB * JStep, n = JB, m=0; n <= JS; n++, J += JStep)
    {
        if ((EBJ = EB - EShift[n]) < 0) EBJ = 0;
        if ((ESJ = ES - EShift[n]) >= NP) ESJ = NP - 1;
        for (i = EBJ; i <= ESJ; i++)
        {
			curVal = Data[n][i] * m_IntensF[n];
            if (curVal < mI) mI = curVal;
            else if (curVal > MI) MI = curVal;
        }
    }
    ISc = 255.999 / (MI - mI) * CF;*/
    for (i=0; i <= JS; ++i) BFJ[i] = -1;
    for (n=0; n < NMarked; ++n)
        if (marked[n].JPos >= JB && marked[n].JPos <= JS
			&& (BFJ[marked[n].JPos] == -1 || abs(marked[n].EPos + EShift[n] - EC) < abs(marked[BFJ[marked[n].JPos]].EPos + EShift[n] - EC)
				|| ((marked[n].LevelState == 0 || marked[n].LevelState == 2 || marked[n].LevelState == 3)
					&&  marked[BFJ[marked[n].JPos]].LevelState != 0 && marked[BFJ[marked[n].JPos]].LevelState != 2 && marked[BFJ[marked[n].JPos]].LevelState != 3)))
			BFJ[marked[n].JPos] = n;
    for (J = JStart + JB * JStep, n = JB, m=0; n <= JS; n++, J += JStep)
	{
		if ((EBJ = EB - EShift[n]) < 0) EBJ = 0;
		if ((ESJ = ES - EShift[n]) >= NP) ESJ = NP - 1;
        for (mI = 1e9, MI = -1e9, i = EBJ; i <= ESJ; i++)
        {
			curVal = Data[n][i];
            if (curVal < mI) mI = curVal;
            else if (curVal > MI) MI = curVal;
        }
        JISc = 255.999 / ((BFJ[n] != -1 ? Data[n][marked[BFJ[n]].EPos] : MI) - mI) * CF;
        for (i = EBJ; i <= ESJ; i++)
		{
			Pixel[n][i] = int((Data[n][i] - mI) * JISc) + intensityOffset;
			if (Pixel[n][i] < 0) Pixel[n][i] = 0;
			else if (Pixel[n][i] > 255) Pixel[n][i] = 255;
		}
		//printf("Pixel[%d][%d]=%d, uchar=%c\n", n, Mind, Pixel[n][Mind], uchar(Pixel[n][Mind]));
	}
	for (m=0; m < NMarked; m++)
		if (marked[m].JPos >= JB && marked[m].JPos <= JS && marked[m].LevelState > 0)
	{
		if ((EBJ = EB - EShift[marked[m].JPos]) < 0) EBJ = 0;
		if ((ESJ = ES - EShift[marked[m].JPos]) >= NP) ESJ = NP - 1;
		if (marked[m].EPos < EBJ || marked[m].EPos > ESJ) continue;
		if (type == Absorption)
		{
			for (m1 = marked[m].EPos - 1; (m1 > EBJ ? 
				Data[marked[m].JPos][m1] < Data[marked[m].JPos][m1 - 1] : false); m1--) ;
			for (m2 = marked[m].EPos + 1; (m2 < ESJ ? 
				Data[marked[m].JPos][m2] < Data[marked[m].JPos][m2 + 1] : false); m2++) ;
		}
		else
		{
			for (m1 = marked[m].EPos - 1; (m1 > EBJ ? 
				Data[marked[m].JPos][m1] > Data[marked[m].JPos][m1 - 1] : false); m1--) ;
			for (m2 = marked[m].EPos + 1; (m2 < ESJ ? 
				Data[marked[m].JPos][m2] > Data[marked[m].JPos][m2 + 1] : false); m2++) ;
		}
		if (m1 < EBJ) m1 = EBJ;
		if (m2 > ESJ) m2 = ESJ;
		if (marked[m].LevelState == 1) 
		{
			for (b = m1; b <= m2; b++) if (Pixel[marked[m].JPos][b] < 256) 
				Pixel[marked[m].JPos][b] += 256;
		}
		else if (marked[m].LevelState >= 2) for (b = m1; b <= m2; b++)
		{
			while (Pixel[marked[m].JPos][b] < 512) Pixel[marked[m].JPos][b] += 256;
			if (marked[m].LevelState >= 4)  Pixel[marked[m].JPos][b] += 256;
		}
	}
	for (n = JB; n <= JS; n++)
	{
		if (EShift[n] > 0)
		{
			for (i = ES, b = ES - EShift[n]; b >= 0 && i >= EB; b--, i--) 
				Pixel[n][i] = Pixel[n][b];
			while (i >= EB) Pixel[n][i--] = 0;
		}
		else if (EShift[n] < 0)
		{
			for (i = EB, b = EB - EShift[n]; b < NP && i <= ES; b++, i++) 
				Pixel[n][i] = Pixel[n][b];
			while (i <= ES) Pixel[n][i++] = 0;
		}
	}
	for (i = EB; i <= ES; i++) 
		for (n = JB, PP = reinterpret_cast<QRgb*>(Pict->scanLine(i)); n <= JS; n++)
	{
		if (Pixel[n][i] < 256) PP[n] = qRgb(Pixel[n][i], Pixel[n][i], Pixel[n][i]);
		else if (Pixel[n][i] < 512) PP[n] = qRgb(0, 0, Pixel[n][i] - 256);
		else if (Pixel[n][i] < 768) PP[n] = qRgb(Pixel[n][i] - 512, 0, 0);
		else PP[n] = qRgb(Pixel[n][i] - 768, Pixel[n][i] - 768, 0);
	}
	Destroy(Pixel, NJ);
	if (Image != 0) delete Image;
	Image = Pict;
	DiagWindow::PSpektrum(P, A, PrintFN);
}

void AddSpectrum::removeLevels()
{
	if (FitD == 0) FitD = State->getFitData();
	int n, N[NPert+1], m;
	TermEnergy TE;
	for (n=0; n <= NPert; n++) N[n] = 0;
	for (n=0; n < NMarked; n++) if (marked[n].LevelState == 3) 
		N[marked[n].StateNum + 1]++;
	int *DelR[NPert+1];
	for (n=0; n <= NPert; n++) if (N[n] > 0) DelR[n] = new int[N[n]];
	for (n=0; n <= NPert; n++) N[n] = 0;
	for (n=0; n < NMarked; n++) if (marked[n].LevelState == 3) 
		DelR[marked[n].StateNum + 1][N[marked[n].StateNum + 1]++] = marked[n].Row;
	for (n=0; (n < NMarked ? marked[n].LevelState < 2 : false); n++) ;
	for (m=n; m < NMarked; m++) 
	{
		if (marked[m].LevelState != markedSaved 
			&& marked[m].LevelState != newAssignment)
		{
			marked[n].JPos = marked[m].JPos;
			marked[n].EPos = marked[m].EPos;
			marked[n].LevelState = marked[m].LevelState;
			marked[n].Row = marked[m].Row;
			marked[n].StateNum = marked[m].StateNum;
			marked[n++].CalcFromProgression = marked[m].CalcFromProgression;
		}
		else if (marked[m].LevelState == 2)
		{
			calcEJ(TE.E, TE.J, marked[m].JPos, marked[m].EPos);
			TE.Iso = Iso;
			TE.v = v;
			TE.err = 0.01;
			TE.dev = 0.0;
			TE.DevR = 0.0;
			TE.ef = (Comp == 0 ? false : true);
			if (marked[m].StateNum == -1 && FitD != 0) 
				FitD->removeMarkedLevel(TE, Source);
			else if (PFitD[marked[m].StateNum] != 0) 
				PFitD[marked[m].StateNum]->removeMarkedLevel(TE, Source);
		}
	}
	NMarked = n;
	for (n=0; n <= NPert; n++) if (N[n] > 0)
	{
		if (n==0 && FitD != 0) FitD->deleteRows(DelR[n], N[n]);
		else if (PFitD[n-1] != 0) PFitD[n-1]->deleteRows(DelR[n], N[n]);
		delete[] DelR[n];
	}
	Paint();
}

void AddSpectrum::setData(double **nData, int nNJ, int nJStart, int nJStep, int nNE, Spektrum *nSource,
                          ElState *nState, int nComp, int nIso, int nv, TermTable *nTT,
                          double nRi, double nRa, double nRes, double *UEnergy, int *vA,
                          int *CA, double *i_IntensityF)
{
	int n;
	Data = nData;
	NMarked = 0;
	NJ = nNJ;
	NP = nNE;
	Source = nSource;
	type = Source->getType();
	Comp = nComp;
	Iso = nIso;
	v = nv;
	TT = (UEnergy == 0 ? nTT : 0);
	NC = (TT != 0 ? TT->getNumComp() : 1);
	TDC = (Comp < NC ? Comp : 0);
    Ri = nRi;
	Ra = nRa;
	Res = nRes;
	if (EShift != 0) delete[] EShift;
	EShift = new int[NJ];
	if (FitD != 0) FitD->clearMarkedLevels();
	for (n=0; n < NPert; n++) if (PFitD[n] != 0) PFitD[n]->clearMarkedLevels();
	JStart = nJStart;
	JStep = nJStep;
    if (m_IntensF != 0) delete[] m_IntensF;
    m_IntensF = i_IntensityF;
	if (UEn != 0) 
	{
        if (vAss != 0) delete[] vAss;
        if (CAss != 0) delete[] CAss;
        vAss = CAss = 0;
        delete[] UEn;
	}
	if (UEnergy == 0)
	{
		double ****TData = TT->getData();
		int J;
		for (n=0, J = JStart, UEn = new double[NJ]; n < NJ; n++, J += JStep) UEn[n] = TData[TDC][Iso][v][J];
	}
	else 
	{
		vAss = vA;
		CAss = CA;
		UEn = UEnergy;
	}
	if (State != nState)
	{
		if (State != 0) disconnect(State, SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
		connect(State = nState, SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
	}
	else setEnergyReference();
	stateChanged();
}

void AddSpectrum::setEnergyReference()
{
	int i, n, m, J, N, MJ, sv = -1, sc = -1;
	bool OFTC;
	double dRes = 1.0 / Res, ****TDB;
	TermTable *TTB = 0;
	Perturbation *Ptn, *PtB;
	ElState *PSt;
	if (OffTT != 0 ? SourceBox->currentText() != OffTT->getName() : true)
	{
		OFTC = true;
		if (SourceBox->currentIndex() == SourceBox->count() - 1)
        {
            if (OffTT == 0) OFTC = false;
            else OffTT = 0;
        }
		else OffTT = State->getTermTable(SourceBox->currentIndex());
		if (PFitD != 0) 
		{
			for (n=0; n < NPert; n++) if (PFitD[n] != 0) 
			{
				disconnect(PFitD[n], SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
				disconnect(PFitD[n], SIGNAL(AssignmentsAccepted(FitData*)), 
						   this, SLOT(AssignmentsAccepted(FitData*)));
				delete[] UCP[n];
			}
			delete[] PFitD;
			delete[] NCP;
			delete[] NvP;
			delete[] NJP;
			delete[] UCP;
			PFitD = 0;
		}
		if (OffTT != 0)
		{
			OTData = OffTT->getData();
			NOC = OffTT->getNumComp();
			NOIso = OffTT->getNumIso();
			NOJ = OffTT->getMaxJ() + 1;
			NOv = OffTT->getMaxv() + 1;
			for (n = NPert = 0, N = OffTT->getNPerturbations(); n<N; n++)
			{
				if ((TTB = OffTT->getPerturbation(n)->Perturber) == 0) continue;
				for(m=0; (m < NPert ? TTB != PTT[m] : false); m++) ;
				if (m == NPert) PTT[NPert++] = TTB;
			}
			for (i=m=0; i < NPert; i++) 
				for (n=0, N = PTT[i]->getNPerturbations(); n<N; n++)
			{
				if ((TTB = PTT[i]->getPerturbation(n)->Perturber) == 0) continue;
				if (TTB != OffTT) for (m=0; (m < NPert ? TTB != PTT[m] : false); m++) ;
				if (m == NPert) PTT[NPert++] = TTB;
			}
			if (NPert > 0)
			{
				PFitD = new FitData*[NPert];
				NCP = new int[NPert];
				NvP = new int[NPert];
				NJP = new int[NPert];
				PTData = new double****[NPert];
				UCP = new bool*[NPert];
				for (n=0; n < NPert; n++)
				{
					PFitD[n] = (PSt = PTT[n]->getElState())->getFitData();
					if (PFitD[n] != 0)
					{
						connect(PFitD[n], SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
						connect(PFitD[n], SIGNAL(AssignmentsAccepted(FitData*)), 
								this, SLOT(AssignmentsAccepted(FitData*)));
					}
					PTData[n] = PTT[n]->getData();
					UCP[n] = new bool[NCP[n] = PTT[n]->getNumComp()];
					NvP[n] = PTT[n]->getMaxv() + 1;
					NJP[n] = PTT[n]->getMaxJ() + 1;
					if (PSt->getLambda() == 1 && NCP[n] == 2)
					{
						UCP[n][0] = true;
						UCP[n][1] = false;
					}
					else for (m=0; m < NCP[n]; m++) UCP[n][m] = true;
				}
			}
			else PFitD = 0;
		}
	}
	else OFTC = false;
    
    if (0 != OffTT && OffTT != TT) N = OffTT->getNPerturbations(Comp, Iso, v) + 3;
    else N = (ResFit != 0 ? ResFit->getNumberOfLocalPerturbations() + 3 : 1);
    if (N==3) N-=2;
    else if (CompBox->count() == 1) CompBox->addItem("least perturbed");
	for (n = CompBox->count() - 1; n>=N; n--) CompBox->removeItem(n);
    for (--n; n < N-2; ++n) CompBox->addItem("part " + QString::number(n));
    if ((N = CompBox->currentIndex()) == 0) for (n=0; n < NJ; n++) EShift[n] = 0;
    else if (0 != OffTT && TT != OffTT)
    {
        if (N == 1) for (n=0, J = JStart; n < NJ; n++, J += JStep)
            EShift[n] = int((UEn[n]
                      - (Iso < NOIso && v < NOv && J < NOJ ? OTData[(Comp < NOC ? Comp : 0)][Iso][v][J] : 0.0)) * dRes);
        else
        {
            if (N>2)
            {
                Ptn = OffTT->getPerturbation(Comp, Iso, v, N-3);
                while (Ptn != 0)
                {
                    m = (TTB = Ptn->Perturber)->getNPerturbations(sc = Ptn->PComp, Iso, sv = Ptn->Pv);
                    for (n=0; n<m; n++)
                    {
                        PtB = TTB->getPerturbation(Ptn->PComp, Iso, Ptn->Pv, n);
                        if (PtB->J >= Ptn->J) break;
                    }
                    if (m>0 && n==m) Ptn = PtB;
                    else if (n>0 && n<m) Ptn = TTB->getPerturbation(Ptn->PComp, Iso, Ptn->Pv, n-1);
                    else Ptn = 0;
                }
            }
            else
            {
                TTB = OffTT;
                sc = Comp;
                sv = v;
            }
            MJ = (Ptn = TTB->getPerturbation(sc, Iso, sv, 0))->J;
            TDB = TTB->getData();
            for (n=0, J = JStart; n < NJ; n++, J += JStep)
            {
                if (J > MJ)
                {
                    if (Ptn == 0) break;
                    TDB = (TTB = Ptn->Perturber)->getData();
                    sc = Ptn->PComp;
                    sv = Ptn->Pv;
                    for (m=0; ((Ptn = TTB->getPerturbation(sc, Iso, sv, m)) != 0 ?
                        Ptn->J < J : false); m++) ;
                    if (Ptn != 0) MJ = Ptn->J;
                    else MJ = TTB->getMaxJ();
                }
                EShift[n] = int((UEn[n] - TDB[sc][Iso][sv][J]) * dRes);
            }
            for ( ; n < NJ; n++) EShift[n] = int(UEn[n] * dRes);
        }
    }
    else
    {
        if (N == 1) for (n=0, J = JStart; n < NJ; n++, J += JStep) EShift[n] = -dRes * ResFit->GetPoint(J);
        else for (n=0, J = JStart; n < NJ; n++, J += JStep) EShift[n] = -dRes * ResFit->GetPoint(J, N-2);
    }
	if (OFTC) stateChanged();
	else Paint();
}

void AddSpectrum::showLevels()
{
	if (NMarked == 0) return;
    for (int n=0; n < NMarked; n++) if (marked[n].LevelState < 1) marked[n].LevelState += 2;
	Paint();
}

void AddSpectrum::stateChanged()
{
	int NL, n, N, m, M, i, NPoints;
	double dRes = 1.0 / Res, E, y;
	TableLine *TL;
	QString TN;
	Point *points = new Point[1000];
	
    /*QFile file("/home/alexst/Documents/Strontium/IntensityComp" + QString::number(vAss != 0 ? vAss[0] : v) + ".dat");
	file.open(QIODevice::WriteOnly);
    QTextStream IntensityCombStream(&file);*/
	
	if (OffTT != 0 ? OffTT->getElState() != State : true)
	{
		SourceBox->blockSignals(true);
		SourceBox->clear();
		for (n=0, N = State->getNumTermTables(); n<N; n++) 
			SourceBox->addItem(State->getTermTableName(n));
		SourceBox->addItem("Potential");
		SourceBox->setCurrentIndex(SourceBox->count() - 1);
		SourceBox->blockSignals(false);
		if (FitD != 0) 
		{
			disconnect(FitD, SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
			disconnect(FitD, SIGNAL(AssignmentsAccepted(FitData*)), 
					   this, SLOT(AssignmentsAccepted(FitData*)));
		}
		if ((FitD = State->getFitData()) != 0)
		{
			connect(FitD, SIGNAL(propertiesChanged()), this, SLOT(stateChanged()));
			connect(FitD, SIGNAL(AssignmentsAccepted(FitData*)), 
					this, SLOT(AssignmentsAccepted(FitData*)));
		}
		for (n=0; n < NJ; n++) EShift[n] = 0;
	}
	else
	{
		SourceBox->blockSignals(true);
		if ((N = State->getNumTermTables()) == SourceBox->count() - 1)
		{
			for (n=0; n<N; n++) if ((TN = State->getTermTableName(n)) != SourceBox->itemText(n))
				SourceBox->setItemText(n, TN);
		}
		else if (N < SourceBox->count())
		{
			for (n=0; (n<N ? 
				State->getTermTableName(n) == (TN = SourceBox->itemText(n)) : false); n++) ;
			SourceBox->removeItem(n);
			if (OffTT != 0 ? TN == OffTT->getName() : true) return;
		}
		else
		{
			for (n=0; (n < SourceBox->count() ? 
				(TN = State->getTermTableName(n)) == SourceBox->itemText(n) : false); n++) ;
			SourceBox->insertItem(n, TN);
		}
		SourceBox->blockSignals(false);
	}
	for (n=0; n < NMarked; n++) 
		if (marked[n].LevelState == -1 || marked[n].LevelState == 3) 
			marked[n].LevelState = 1;
	for (n=0; (n < NMarked ? marked[n].LevelState != 1 : false); n++) ;
	for (m=n+1; m < NMarked; m++) if (marked[m].LevelState != 1) 
		marked[n++] = marked[m];
	NMarked = n;
	for(n=0; n < NMarked; n++)
	{
		points[n].x = marked[n].JPos * JStep + JStart;
		E = (double(marked[n].EPos) - 0.5 * NP) * Res;
		ParabInterpol(E - Res, Data[marked[n].JPos][marked[n].EPos - 1], E,
						  Data[marked[n].JPos][marked[n].EPos], E + Res,
						  Data[marked[n].JPos][marked[n].EPos + 1], points[n].y, y);
		points[n].unc = 0.0201;
		points[n].sig = 4975.12437811;
		points[n].MarkerNum = n;
	}
	for (i=-1, NPoints = NMarked; i < NPert; i++)
	{
		if (i==-1 && FitD != 0) FitD->getData(TL, NL);
		else if (i>=0 ? PFitD[i] != 0 : false) PFitD[i]->getData(TL, NL);
		else continue;
		for (n=N=0; n < NL && NMarked < 1000; n++) 
			if (TL[n].Iso == Iso && ((TL[n].Js != TL[n].Jss && Comp == 0) 
				|| (TL[n].Js == TL[n].Jss && Comp == 1))
				&& (m = (TL[n].Jss - JStart) / JStep) < NJ) 
		{
			M = int(0.5 * NP + (TL[n].WN - UEn[m]) * dRes);
			if (M < 0 || M >= NP) continue;
			M = getPeakTop(m, M);
			marked[NMarked].JPos = m;
			marked[NMarked].EPos = M;
			marked[NMarked].Row = TL[n].Row;
			marked[NMarked].LevelState = 1;
			marked[NMarked].CalcFromProgression = (TL[n].LTab != 0);
			marked[NMarked].CorrectVibLevel = (TL[n].vss == (vAss != 0 ? vAss[TL[n].Jss] : v));
			if (marked[NMarked].CorrectVibLevel)
			{
				points[NPoints].x = TL[n].Jss;
				points[NPoints].y = TL[n].WN - UEn[m];
				points[NPoints].unc = TL[n].err;
                points[NPoints].obs = TL[n].WN;
				points[NPoints].sig = 1.0 / TL[n].err;
				points[NPoints++].MarkerNum = NMarked;
			}
			else marked[NMarked].SplineNum = -1;
			marked[NMarked++].StateNum = i;
            //IntensityCombStream << QString("%1\t%2\t%3\t%4\n").arg(TL[n].Jss).arg(M).arg(Data[m][M], 0, 'e', 6).arg(m_IntensF[m], 0, 'e', 6);
		}
		delete[] TL;
	}
	if (NPoints > 0)
	{
        ResidualFit* newResidualFit = FitD->getResidualFit(State, Iso, v, Comp);
        if (newResidualFit != ResFit)
        {
            if (newResidualFit != 0)
            {
                if (ResFit != 0 && m_ownResidualFit)
                {
                    delete ResFit;
                    m_ownResidualFit = false;
                }
                Molecule *mol = State->getMolecule();
                IsoTab* isoTab = mol->getIso();
                ResFit = newResidualFit;
                ResFit->setFitData(NPoints, points, isoTab->relRedMass[Iso], State->getOmega(), mol->getJStep(Iso), State->getBe());
                setEnergyReference();
                delete isoTab;
            }
            else
            {
                if (ResFit == 0 || !m_ownResidualFit)
                {
                    ResFit = new ResidualFit(State);
                    m_ownResidualFit = true;
                }
                ResFit->FitTo(NPoints, points);
                for (n=0; n < NPoints; n++) 
                    marked[points[n].MarkerNum].SplineNum = points[n].SplineNum;
                drawDAddSpectrum(dRes, true);
            }
        }
	}
	//printf("NMarked = %d\n", NMarked);
	Paint();
}

int AddSpectrum::getPeakTop(int m, int M)
{
	if (type == Absorption) 
	{
		while (M>0 ? Data[m][M-1] < Data[m][M] : false) M--;
		while (M < NP - 1 ? Data[m][M+1] < Data[m][M] : false) M++;
	}
	else
	{
		while (M>0 ? Data[m][M-1] > Data[m][M] : false) M--;
		while (M < NP - 1 ? Data[m][M+1] > Data[m][M] : false) M++;
	}
	return M;
}

void AddSpectrum::drawDAddSpectrum(double dRes, bool newState, int SplineNum, int JMin, int JMax)
{
	double **DAddData = Create(NP, 2), E, F1, F2;
	int J, d, n, m, M, i, N = ResFit->getNumSplines();
	Spline *spline;
	for (n=0, E = -0.5 * double(NP) * Res; n < NP; n++, E += Res) 
	{
		DAddData[n][0] = E;
		DAddData[n][1] = 0.0;
	}
	if (SplineNum >= 0 && SplineNum < N) N = SplineNum + 1;
	else SplineNum = 0;
	for (n = SplineNum; n<N; n++)
	{
		spline = ResFit->getSpline(n);
		for (J = ((J = int(spline->getx0())) < JMin ? JMin : J), M = ((M = int(spline->getxN())) > JMax && JMax >= 0 ? JMax : M); J <= M; J += JStep)
		{
			m = (J - JStart) / JStep; 
			E = spline->gety(J) * dRes;
			d = int(E);
			F2 = E - double(d);
			F1 = 1.0 - F2;
			if (d>=0) for (i=0; d+1 < NP; i++, d++) 
				DAddData[i][1] += F1 * Data[m][d] + F2 * Data[m][d+1];
			else for (i=-d; i < NP; i++)
				DAddData[i][1] += F1 * Data[m][i+d] + F2 * Data[m][i+d+1];
		}
	}
	if (Source->getType() == Absorption) for (n=0; n < NP; n++) DAddData[n][1] *= -1;
	if (newState)
	{
		if (DAddSpect == 0)
		{
			DAddSpect = new DoubleAddSpectrum(this, MW);
			MW->showMDIChild(DAddSpect);
		}
		DAddSpect->setData(DAddData, NP);
	}
	else DAddSpect->addData(DAddData, NP);
	Destroy(DAddData, NP);
}

void AddSpectrum::maximizeDAdd(double Diff)
{
	if (ResFit == 0) return;
	double EStart = -0.5 * Res * double(NP), dRes = 1.0 / Res;
	DAddMaximize(ResFit, Data, NJ, JStart, JStep, NP,  EStart, dRes, Diff, Source->getType() == Absorption);
	drawDAddSpectrum(dRes, false);
}

int AddSpectrum::AutoUpdateLevels(bool fullyAutomatic)
{
	int j, n, N =  (ResFit != 0 ? ResFit->getNumSplines() : 0), M, m, Result = autoUpdateSuccess;
	if (NMarked == 0) m = M = 0;
	else 
	{
		for (n=0; n < NMarked && !marked[n].CorrectVibLevel; n++) ;
		if (marked[n].CorrectVibLevel) 
		{
			for (m = M = marked[n++].JPos; n < NMarked; n++) if (marked[n].CorrectVibLevel)
			{
				if (marked[n].JPos < m) m = marked[n].JPos;
				else if (marked[n].JPos > M) M = marked[n].JPos;
			}
		}
		else m=M=0;
	}
	if (M-m <= JStep)
	{
		if (!fullyAutomatic)
			QMessageBox::information(this, "MolSpektanalysis", 
								 "For this vibrational state the available amount of \
								 initial data is insufficient for the current \
								 implementation of this function!");
		return autoUpdateInsufficientData;
	}
	double Ec, E, y;
	double Unc = FitD->getUncertaintyOfvibLevel((vAss != 0 ? vAss[0] : v), Iso, Source);
	double dRes = 1.0 / Res, diff = DAddSpect->getDiff();
	if (isnan(diff) || isinf(diff))
	{
		if (!fullyAutomatic) QMessageBox::information(this, "MolSpektanalysis", "Invalid double add spectrum!");
		return autoUpdateCompleteFailure;
	}
	Spline *sp;
	NMBAUL = NMarked;
	maximizeDAdd(diff);
	for (n=M=0; n<N; n++)
	{
		sp = ResFit->getSpline(n);
		m = (sp->getxN() - sp->getx0()) / JStep + 1;
		if (m>M) M=m;
	}
	int **Peaks = CreateInt(N, M), **JPos = CreateInt(N, M);
	int MJ[N];
	bool OK;
	double Eo, SNR;
	TermEnergy TE;
	for (n=0; n<N; n++)
	{
		for (sp = ResFit->getSpline(n), j = int(sp->getx0()), m=0; j <= sp->getxN(); 
			 j += JStep, m++)
		{
			if (m==0) JPos[n][0] = (j - JStart) / JStep;
			else JPos[n][m] = JPos[n][m-1] + 1;
			Peaks[n][m] = int(0.5 * NP + (Ec = sp->gety(j)) * dRes);
			if (Peaks[n][m] > 0 && Peaks[n][m] < NP - 1) 
				Peaks[n][m] = getPeakTop(JPos[n][m], Peaks[n][m]);
			else 
			{
				Peaks[n][m] = -1;
				continue;
			}
			E = (double(Peaks[n][m]) - 0.5 * NP) * Res;
			ParabInterpol(E - Res, Data[JPos[n][m]][Peaks[n][m] - 1], E,
						  Data[JPos[n][m]][Peaks[n][m]], E + Res,
						  Data[JPos[n][m]][Peaks[n][m] + 1], Eo, y);
			if (fabs(Ec - Eo) > 2.0 * Unc) Peaks[n][m] = -1;
		}
		MJ[n] = m-1;
	}
	for (n=0, M = NMarked; n<M; n++) if (marked[n].CorrectVibLevel)
	{
		if (marked[n].LevelState <= 0) marked[n].LevelState += 2;
		else if (marked[n].LevelState == 3) marked[n].LevelState = 1;
		for (m=0, OK = false; m<N; m++) 
			if (marked[n].JPos >= JPos[m][0] && marked[n].JPos <= JPos[m][MJ[m]])
		{
			j = marked[n].JPos - JPos[m][0];
			if (Peaks[m][j] > 0)
			{
				if (marked[n].EPos == Peaks[m][j])
				{
					Peaks[m][j] = -2;
					OK = true;
					break;
				}
			}
			else if (Peaks[m][j] == -1)
			{
				E = (double(marked[n].EPos) - 0.5 * NP) * Res;
				ParabInterpol(E - Res, Data[marked[n].JPos][marked[n].EPos - 1], E,
							  Data[marked[n].JPos][marked[n].EPos], E + Res,
							  Data[marked[n].JPos][marked[n].EPos + 1], Ec, y);
				if (fabs(Ec - ResFit->getSpline(m)->gety(j)) < 4.0 * Unc) OK = true;
			}
		}
		if (!OK) marked[n].LevelState += 3;
	}
	int NP[N], NRem[N], k, l, s;
	bool Failure[N];
	for (n=0; n<N; n++) NP[n] = NRem[n] = 0;
	/*for (n=0; n < NMarked; n++) if (marked[n].CalcFromProgression)
	{
		NP[marked[n].SplineNum]++;
		if (marked[n].LevelState == savedToDelete 
			|| marked[n].LevelState == newToDelete) NRem[marked[n].SplineNum]++;
	}
	for (n=0; n<N; n++)
	{
		if (NP[n] == 0) Failure[n] = true;
		else Failure[n] = false;
	}*/
	for (n=0; n < NMarked; n++) if (marked[n].CorrectVibLevel) //if (Failure[marked[n].SplineNum])
	{
		if (marked[n].LevelState == savedAssignment) NP[marked[n].SplineNum]++;
		if (marked[n].LevelState == savedToDelete) NRem[marked[n].SplineNum]++;
	}
	for (n=0; n<N; n++) 
	{
		Failure[n] = (NRem[n] != 0 && double(NP[n]) / NRem[n] <= 3.0);
		NP[n] = 0;
	}
	for (n=0; n<N; n++) for (m=0; m <= MJ[n]; m++) if (Peaks[n][m] > 0)
	{
		fillTEStruct(TE, j, JPos[n][m], Peaks[n][m]);
		addLevel(TE, j, JPos[n][m], Peaks[n][m]);
		marked[NMarked - 1].SplineNum = n;
		marked[NMarked - 1].CorrectVibLevel = true;
	}
	Destroy(Peaks, N);
	Destroy(JPos, N);
	if (fullyAutomatic)
	{
		MarkedPeak MB;
		for (MB.JPos = 1; MB.JPos != -1; ) for (n=1, MB.JPos = -1; n < NMarked; n++) 
			if (marked[n].SplineNum < marked[n-1].SplineNum || (marked[n].SplineNum == marked[n-1].SplineNum && marked[n].JPos < marked[n-1].JPos))
		{
			MB = marked[n-1];
			marked[n-1] = marked[n];
			marked[n] = MB;
		}
		for (m=0; m < NMarked && marked[m].SplineNum == -1; m++) ;
		for (n=l=0; n<N; n++)
		{ 
			for (j=k=0, l=m; m < NMarked && marked[m].SplineNum == n; m++)
			{
				if (marked[m].LevelState == newAssignment) j++;
				if (marked[m].LevelState == savedAssignment) k++;
			}
			if (j>k || Failure[n])
			{
				drawDAddSpectrum(dRes, true, n, marked[l].JPos, marked[m-1].JPos);
				SNR = DAddSpect->getSNR();
				if (SNR < 1.0 || (Failure[n] && SNR < sqrt(double(j+k)*0.1)))
				{
					Failure[n] = true;
					for (j=l; j<m; j++) 
					{
						if (marked[j].LevelState == newAssignment) marked[j].LevelState = newToDelete;
						else if (marked[j].LevelState == savedToDelete) marked[j].LevelState = savedAssignment;
					}
				}
				else Failure[n] = false;
			}
			if (!Failure[n]) for (j=l, k=0; j<=m; j++)
			{
				if (k==0) s=j;
				if (j<m && marked[j].LevelState == newAssignment) k++;
				if (j==m || marked[j].LevelState == savedAssignment)
				{
					if (k >= 10)
					{
						drawDAddSpectrum(dRes, true, n, marked[s].JPos, marked[j-1].JPos);
						if (DAddSpect->getSNR() < 1.0) for (k=s; k<j; k++) 
						{
							if (marked[k].LevelState == newAssignment) marked[k].LevelState = newToDelete;
							else if (marked[k].LevelState == savedToDelete) marked[k].LevelState = savedAssignment;
						}
					}
					k=0;
				}
			}
		}
	}
	for (n=0; n < NMarked; n++) 
		if ((marked[n].CalcFromProgression || (marked[n].CorrectVibLevel && Failure[marked[n].SplineNum])) && marked[n].LevelState == savedToDelete) 
		marked[n].LevelState = savedAssignment;
	for (n=0; n<N && Failure[n]; n++) ;
	if (n==0)
	{
		for (n=0; n<N && !Failure[n]; n++) ;
		if (n==N) Result = autoUpdateSuccess;
		else Result = autoUpdatePartialFailure;
	}
	else if (n==N) Result = autoUpdateCompleteFailure;
	else Result = autoUpdatePartialFailure;
	Paint();
	return Result;
}

void AddSpectrum::acceptAutoUpdateResult()
{
	int n;
	for (n=0; n < NMarked; n++)
	{
		if (marked[n].LevelState == newAssignment) marked[n].LevelState = hiddenNew;
		else if (marked[n].LevelState == newToDelete) marked[n].LevelState = newAssignment;
		else if (marked[n].LevelState == savedToDelete) marked[n].LevelState = markedSaved;
	}
	removeLevels();
	updateLevels();
	FitD->writeData();
	for (n=0; n < NPert; n++) PFitD[n]->writeData();
	NMBAUL = -1;
	Paint();
}

void AddSpectrum::rejectAutoUpdateResult()
{
	int n;
	for (n=0; n < NMBAUL; n++) if (marked[n].LevelState == newAssignment)
		marked[n].LevelState = newToDelete;
	removeLevels();
	for (n=0; n < NMarked; n++)
	{
		if (marked[n].LevelState == newToDelete) marked[n].LevelState = newAssignment;
		else if (marked[n].LevelState == savedToDelete) marked[n].LevelState = savedAssignment;
	}
	NMBAUL = -1;
	Paint();
}

void AddSpectrum::updateLevels()
{
	int c, n, J, N[NPert+1], **r = CreateInt(NPert + 1, NMarked), v[NMarked], C[NMarked];
	double **E = Create(NPert + 1, NMarked);
	if (FitD == 0) FitD = State->getFitData();
	for (n=0; n <= NPert; n++) N[n] = 0;
	for (n=0; (n < NMarked ? marked[n].LevelState < 3 : false); n++) ;
	if (n < NMarked) 
	{
		for (n=0; n < NMarked; n++) if (marked[n].LevelState == 3)
		{
			c = marked[n].StateNum + 1;
			calcEJ(E[c][N[c]], J, marked[n].JPos, marked[n].EPos);
			if (vAss != 0)
			{
				v[N[c]] = vAss[marked[n].JPos];
				C[N[c]] = -1 - CAss[marked[n].JPos];
			}
			r[c][N[c]++] = marked[n].Row;
		}
	}
	else for (n=0; n < NMarked; n++) if (marked[n].LevelState == 1)
	{
		c = marked[n].StateNum + 1;
		calcEJ(E[c][N[c]], J, marked[n].JPos, marked[n].EPos);
		if (vAss != 0)
		{
			v[N[c]] = vAss[marked[n].JPos];
			C[N[c]] = -1 - CAss[marked[n].JPos];
		}
		r[c][N[c]++] = marked[n].Row;
	}
	if (N[0] > 0) 
	{
		if (NPert > 0 || vAss == 0) FitD->updateEnergy(N[0], r[0], E[0]);
		else FitD->updateLevels(N[0], r[0], E[0], C, v);
	}
	for (n=1; n <= NPert; n++) if (N[n] > 0) PFitD[n-1]->updateEnergy(N[n], r[n], E[n]);
	Destroy(r, NPert + 1);
	Destroy(E, NPert + 1);
}
