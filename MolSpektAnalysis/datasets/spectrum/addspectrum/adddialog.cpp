//
// C++ Implementation: adddialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#include "adddialog.h"
#include "MainWindow.h"
#include "Spektrum.h"
#include "elstate.h"
#include "molecule.h"
#include "constants.h"
#include "utils.h"
#include "termtable.h"
#include "fitdata.h"
#include "linetable.h"
#include "potential.h"
#include "naturalspline.h"
#include "isotab.h"
#include "CoupledSineWaveFunc.h"
#include "addspektrum.h"

#include <stdio.h>
#include <math.h>

#include <QPushButton>
#include <QRadioButton>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QGridLayout>
#include <QMessageBox>
#include <QCheckBox>
#include <QMenu>
#include <QToolButton>
#include <QTextEdit>
#include <QTextStream>


AddDialog::AddDialog(MainWindow *mw, Spektrum *nSpekt) : QWidget()
{
	MW = mw;
	Image = 0;
	TData = 0;
	IData = 0;
	NJ = NC = 0;
	LT = TT = 0;
	Spekt = nSpekt;
	Marked = false;
	LogWindow = 0;
	int n, M = MW->getNumMolecules();
	setAttribute(Qt::WA_DeleteOnClose);
	if (mw == 0)
	{
		printf("AddDialog error: mw==0, closing!!");
		close();
	}
	setWindowTitle("Add shifted versions of a spectrum");
	QGridLayout *Layout = new QGridLayout(this);
	QLabel *L3 = new QLabel("Spectrum:", this);
	Layout->addWidget(L3, 0, 0);
	SpektBox = new QComboBox(this);
	SpektBox->setEditable(false);
	Layout->addWidget(SpektBox, 0, 1);
	QLabel *L4 = new QLabel("Molecule:", this);
	Layout->addWidget(L4, 0, 2);
	MolBox = new QComboBox(this);
	MolBox->setEditable(false);
	Layout->addWidget(MolBox, 0, 3);
	QLabel *L5 = new QLabel("Lower state:", this);
	Layout->addWidget(L5, 0, 4);
	LSBox = new QComboBox(this);
	Layout->addWidget(LSBox, 0, 5);
	LSBox->setEditable(false);
	QLabel *L6 = new QLabel("Upper state:", this);
	Layout->addWidget(L6, 0, 6);
	USBox = new QComboBox(this);
	USBox->setEditable(false);
	Layout->addWidget(USBox, 0, 7);
	QGridLayout *GL2 = new QGridLayout;
	AddByBand = new QRadioButton("Add by band", this);
	GL2->addWidget(AddByBand, 0, 0);
	AddByProg = new QRadioButton("Add by progression", this);
	GL2->addWidget(AddByProg, 0, 1);
	AddByProgFCF = new QRadioButton("Add by prog FCF", this);
	GL2->addWidget(AddByProgFCF, 0, 3);
	AutoAnalyse = new QCheckBox("Search all levels", this);
	GL2->addWidget(AutoAnalyse, 0, 4);
	connect(AutoAnalyse, SIGNAL(stateChanged(int)), this, SLOT(MethodChanged()));
	DrawImage = new QCheckBox("Draw J dep. Image", this);
	DrawImage->setEnabled(false);
	GL2->addWidget(DrawImage, 0, 5);
	connect(DrawImage, SIGNAL(stateChanged(int)), this, SLOT(MethodChanged()));
	UpdateAllVs = new QCheckBox("Auto update all levels", this);
	UpdateAllVs->setEnabled(false);
	connect(UpdateAllVs, SIGNAL(stateChanged(int)), this, SLOT(MethodChanged()));
	GL2->addWidget(UpdateAllVs, 0, 6);
	Layout->addLayout(GL2, 1, 0, 1, 8);
	QLabel *L7 = new QLabel("Isotopologue:", this);
	Layout->addWidget(L7, 2, 0);
	IsoBox = new QComboBox(this);
	IsoBox->setEditable(false);
	Layout->addWidget(IsoBox, 2, 1);
	L1 = new QLabel("v':", this);
	Layout->addWidget(L1, 2, 2);
	L1E = new QLineEdit("0", this);
	Layout->addWidget(L1E, 2, 3);
	L2 = new QLabel("J':", this);
	Layout->addWidget(L2, 2, 4);
	L2E = new QLineEdit("15", this);
	Layout->addWidget(L2E, 2, 5);
	QLabel *L8 = new QLabel("Prog type:", this);
	Layout->addWidget(L8, 2, 6);
	CompBox = new QComboBox(this);
	CompBox->addItem("Q lines");
	CompBox->addItem("P/R doublets");
	CompBox->setCurrentIndex(1);
	Layout->addWidget(CompBox, 2, 7);
	QLabel *L9 = new QLabel("Ref. spectrum:", this);
	Layout->addWidget(L9, 3, 0);
	RSpektBox = new QComboBox(this);
	Layout->addWidget(RSpektBox, 3, 1);
	QLabel *L10 = new QLabel("Ref. iso:", this);
	Layout->addWidget(L10, 3, 2);
	RIsoBox = new QComboBox(this);
	Layout->addWidget(RIsoBox, 3, 3);
	QLabel *L11 = new QLabel("Ref. J':", this);
	Layout->addWidget(L11, 3, 4);
	RJsE = new QLineEdit("15", this);
	Layout->addWidget(RJsE, 3, 5);
	QLabel *L12 = new QLabel("Ref. prog.:", this);
	Layout->addWidget(L12, 3, 6);
	RCompBox = new QComboBox(this);
	RCompBox->addItem("Q lines");
	RCompBox->addItem("P/R doublets");
	RCompBox->setCurrentIndex(1);
	Layout->addWidget(RCompBox, 3, 7);
	QGridLayout *TL = new QGridLayout;
	QLabel *L13 = new QLabel("Temperature:", this);
	TL->addWidget(L13, 0, 0);
	TempE = new QLineEdit("950", this);
	TL->addWidget(TempE, 0, 1);
	TUnit = new QComboBox(this);
	TUnit->addItem(QString("°C").right(2));
	TUnit->addItem("K");
	TUnit->setEditable(false);
	TL->addWidget(TUnit, 0, 2);
	Layout->addLayout(TL, 4, 0, 1, 2);
	L14 = new QLabel("R min [cm^-1]:", this);
	Layout->addWidget(L14, 4, 2);
	RiE = new QLineEdit("18000", this);
	Layout->addWidget(RiE, 4, 3);
	QLabel *L15 = new QLabel("R max [cm^-1]:", this);
	Layout->addWidget(L15, 4, 4);
	RaE = new QLineEdit("18020", this);
	Layout->addWidget(RaE, 4, 5);
	QLabel *L16 = new QLabel("Resolution [cm^-1]:", this);
	Layout->addWidget(L16, 4, 6);
	ResE = new QLineEdit("0.01", this);
	Layout->addWidget(ResE, 4, 7);
	Layout->setRowMinimumHeight(20, 5);
	QGridLayout *BL = new QGridLayout;
	ApplyB = new QPushButton("Apply", this);
	connect(ApplyB, SIGNAL(pressed()), this, SLOT(Apply()));
	BL->addWidget(ApplyB, 0, 0);
	CloseB = new QPushButton("Close", this);
	connect(CloseB, SIGNAL(pressed()), this, SLOT(close()));
	BL->addWidget(CloseB, 0, 1);
	Layout->addLayout(BL, 6, 0, 1, 8);
	connect(MolBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(MolChanged(QString)));
	for (n=0; n<M; n++) MolBox->addItem(MW->getMolecule(n)->getName());
	connect(AddByBand, SIGNAL(toggled(bool)), this, SLOT(MethodChanged()));
	connect(AddByProg, SIGNAL(toggled(bool)), this, SLOT(MethodChanged()));
	AddByBand->setChecked(true);
	connect(MW, SIGNAL(SpectrumChanged(Spektrum*)), this, SLOT(UpdateSpektList()));
	connect(MW, SIGNAL(SpektrumDeleted(int)), this, SLOT(UpdateSpektList()));
	UpdateSpektList();
}

AddDialog::~AddDialog()
{
	
}

void AddDialog::enableOptimize(bool)
{
	/*if ((Marked = enable) && DrawImage->isChecked()) Optimize->setEnabled(true);
	else Optimize->setEnabled(false);*/
}

void AddDialog::UpdateSpektList()
{
	int n, N = MW->getNumSpectra();
	QString SN, RSN = RSpektBox->currentText();
	Spektrum *SB;
	SpektBox->clear();
	RSpektBox->clear();
	for (n=0; n<N; n++)
	{
		SN = (SB = MW->getSpectrum(n))->getName();
		SpektBox->addItem(SN);
		RSpektBox->addItem(SN);
		if (Spekt == SB)
		{
			SpektBox->setCurrentIndex(n);
			if (RSN.isEmpty()) RSpektBox->setCurrentIndex(n);
		}
		if (SN == RSN) RSpektBox->setCurrentIndex(n);
	}
}

void AddDialog::MolChanged(QString Name)
{
	Molecule *Mol = MW->getMolecule(Name);
	int n, N = Mol->getNumStates();
	QString SN;
	LSBox->clear();
	USBox->clear();
	IsoBox->clear();
	RIsoBox->clear();
	for (n=0; n<N; n++)
	{
		SN = Mol->getState(n);
		LSBox->addItem(SN);
		USBox->addItem(SN);
	}
	IsoTab *Iso = Mol->getIso();
	for (n=0, NumIso = Iso->numIso; n < NumIso; n++)
	{
		IsoBox->addItem(SN = Iso->getIsoName(n));
		RIsoBox->addItem(SN);
	}
	UpdateAllVs->setEnabled(false);
}

void AddDialog::MethodChanged()
{
	bool TDA = (MW->getMolecule(MolBox->currentText())->getStateP(USBox->currentIndex())->getTermTable() == 0 ? 
				false : true);
	if (AutoAnalyse->isChecked())
	{
		DrawImage->setEnabled(false);
		UpdateAllVs->setEnabled(false);
	}
	else if (DrawImage->isChecked())
	{
		AutoAnalyse->setEnabled(false);
		UpdateAllVs->setEnabled(true);
	}
	else
	{
		DrawImage->setEnabled(true);
		if (AddByProgFCF->isChecked()) AutoAnalyse->setEnabled(true);
		else AutoAnalyse->setEnabled(false);
		UpdateAllVs->setEnabled(false);
	}
	if (UpdateAllVs->isChecked() && UpdateAllVs->isEnabled())
	{
		if (CompBox->count() == 2) CompBox->addItem("All components");
		AddByBand->setEnabled(false);
		AddByProg->setEnabled(false);
		AddByProgFCF->setEnabled(true);
		RCompBox->setEnabled(false);
		RIsoBox->setEnabled(false);
		RJsE->setEnabled(false);
		RSpektBox->setEnabled(false);
		TempE->setEnabled(true);
		TUnit->setEnabled(true);
		USBox->setEnabled(true);
		LSBox->setEnabled(true);
		MolBox->setEnabled(true);
		if (IsoBox->count() == NumIso) IsoBox->addItem("All iso");
		IsoBox->setCurrentIndex(NumIso);
		IsoBox->setEnabled(true);
		L1->setText("v min:");
		L2->setText("v max:");
		L1E->setEnabled(true);
		L2E->setEnabled(true);
	}
	else
	{
		if (IsoBox->count() > NumIso) IsoBox->removeItem(NumIso);
		if (CompBox->count() == 3) CompBox->removeItem(2);
		MolBox->setEnabled(true);
		LSBox->setEnabled(true);
		IsoBox->setEnabled(true);
		AddByBand->setEnabled(true);
		AddByProg->setEnabled(true);
		AddByProgFCF->setEnabled(true);
		if (AddByBand->isChecked())
		{
			L1->setText("v':");
			L2->setText("v'':");
			L1E->setEnabled(true);
			L2E->setEnabled(true);
			USBox->setEnabled(true);
			TempE->setEnabled(true);
			TUnit->setEnabled(true);
			RCompBox->setEnabled(false);
			RIsoBox->setEnabled(false);
			RJsE->setEnabled(false);
			RSpektBox->setEnabled(false);
			if (RaE->isEnabled() && TDA)
			{ 
				RiE->setText(QString::number((RaE->text().toDouble() - RiE->text().toDouble()) / 2.0, 'g', 8));
				L14->setText("Search radius:");
				RaE->setEnabled(false);
			}
			else if (!RaE->isEnabled() && !TDA)
			{
				RiE->setText(QString::number(RaE->text().toDouble() - 2.0 * RiE->text().toDouble(), 'g', 8));
				L14->setText("R min [cm^-1]:");
				RaE->setEnabled(true);
			}
		}
		else if (AddByProg->isChecked())
		{
			TempE->setEnabled(false);
			TUnit->setEnabled(false);
			RCompBox->setEnabled(true);
			RIsoBox->setEnabled(true);
			RJsE->setEnabled(true);
			RSpektBox->setEnabled(true);
			L2E->setEnabled(false);
			if (DrawImage->isChecked())
			{
				L1->setText("v':");
				L1E->setEnabled(true);
				USBox->setEnabled(true);
				if (RaE->isEnabled())
				{ 
					RiE->setText(QString::number((RaE->text().toDouble() - RiE->text().toDouble()) / 2.0, 'g', 8));
					L14->setText("Search radius:");
					RaE->setEnabled(false);
				}
			}
			else
			{
				L1->setText("J':");
				L1E->setEnabled(true);
				USBox->setEnabled(false);
				if (!RaE->isEnabled())
				{
					RiE->setText(QString::number(RaE->text().toDouble() - 2.0 * RiE->text().toDouble(), 'g', 8));
					L14->setText("R min [cm^-1]:");
					RaE->setEnabled(true);
				}
			}
		}
		else if (AddByProgFCF->isChecked())
		{
			if (AutoAnalyse->isChecked())
			{
				L1E->setEnabled(false);
				L2E->setEnabled(false);
			}
			else
			{
				L1->setText("v':");
				L1E->setEnabled(true);
				if (DrawImage->isChecked()) L2E->setEnabled(false);
				else
				{
					L2->setText("J':");
					L2E->setEnabled(true);
				}
			}
			USBox->setEnabled(true);
			TempE->setEnabled(true);
			TUnit->setEnabled(true);
			RCompBox->setEnabled(false);
			RIsoBox->setEnabled(false);
			RJsE->setEnabled(false);
			RSpektBox->setEnabled(false);
			if (RaE->isEnabled())
			{ 
				RiE->setText(QString::number((RaE->text().toDouble() - RiE->text().toDouble()) / 2.0, 'g', 8));
				L14->setText("Search radius:");
				RaE->setEnabled(false);
			}
		}
	}
}

void AddDialog::Apply()
{
	double T = TempE->text().toDouble();
	double Res = ResE->text().toDouble();
	double Ri = RiE->text().toDouble();
	double Ra = RaE->text().toDouble();
    int n, N = MW->getNumSpectra(), v, I = -1, NumWFPoints = NumPoints;
	int NP = int((RaE->isEnabled() ? Ra - Ri : 2.0 * Ri) / Res) + 1;
	double **Result = Create(NP, 2);
	bool S = false;
	QString SN = SpektBox->currentText(), RSN = RSpektBox->currentText(), BSN;
	Spektrum *RSpekt = 0, *SB;
	Molecule *Mol = MW->getMolecule(MolBox->currentText());
	//printf("AddDialog::Apply: Result=%d, NP=%d\n", Result, NP);
	if (Mol == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: The selected molecule doesn't exist anymore!", QMessageBox::Ok);
		return;
	}
	if (TUnit->currentIndex() == 0) T += C_K_C;
	for (Spekt = 0, n=0; n<N; n++) 
	{
		BSN = (SB = MW->getSpectrum(n))->getName();
		if (BSN == SN) Spekt = SB;
		if (BSN == RSN) RSpekt = SB;
	}
	if (Spekt == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: The selected spectrum doesn't exist anymore!", QMessageBox::Ok);
		return;
	}
	ElState *UState = Mol->getStateP(USBox->currentIndex());
	ElState *LState = Mol->getStateP(LSBox->currentIndex());
	if (UState->getName() != USBox->currentText() || LState->getName() != LSBox->currentText())
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: A selected state changed its name!", QMessageBox::Ok);
		return;
	}
	if (DrawImage->isChecked())
	{
        double ****LTData = (LT != 0 ? LT->getData() : 0), *WF = 0, *Off = 0, Rd = Ri, ***UWF = 0, *UEn = 0, TStr[10], *IntensF = 0;
		int MJ, i, j, k, l, m, c, J, Mv, NOff = 0, RIso = -1, RComp = -1, RJs = -1, vss, Comp, Iso, *JA = 0, *vA = 0, *CA = 0, WFC[10], NWFC = 0;
		int vmin = L1E->text().toInt(), vmax = L2E->text().toInt(), iso = IsoBox->currentIndex();
		if (iso == NumIso) iso = -1;
		//int CompMin, CompMax, IsoMin, IsoMax, vMin, vMax;
		//bool Opt = Optimize->isEnabled() && Optimize->isChecked();
		CoupledSineWaveFunc *SWF = (UState->getPotential() != 0 ? UState->getPotential()->getCoupledSineWaveFunc() : 0);
		if (AddByProg->isChecked())
		{
			if (RSpekt == 0)
			{
				QMessageBox::information(this, "MolSpektanalysis", 
						"Error: The selected reference spectrum doesn't exist anymore!", QMessageBox::Ok);
				return;
			}
			RIso = RIsoBox->currentIndex();
			RComp = RCompBox->currentIndex();
			RJs = RJsE->text().toInt();
		}
		/*if (Opt)
		{
			return;
			//Muss überarbeitet werden, bevor funktionsfähig
			int NMarked, **marked;
			Image->getMarked(NMarked, marked);
			c = (Comp < NC ? Comp : 0);
			for (i=1, j = marked[0][0]; i < NMarked; i++) if (marked[i][0] > j) j = marked[i][0];
			for (Nv = 0, J = JStart + j * JStep, Mv = LT->getMaxv(); (Nv <= Mv ? LTData[0][Iso][Nv][J] != 0.0 : false); 
				 Nv++) ;
			double E, **Offset = Create(NMarked, 2 * Nv), TE[NMarked], y;
			WF = new double[2 * Nv];
			for (i=0; i < NMarked; i++)
			{
				J = JStart + marked[i][0] * JStep;
				E = TData[c][Iso][v][J] - Ri + double(marked[i][1]) * Res;
				if (marked[i][1] > 0 && marked[i][1] < NP - 1) 
					ParabInterpol(E - Res, IData[marked[i][0]][marked[i][1] - 1], E, IData[marked[i][0]][marked[i][1]], 
								  E + Res, IData[marked[i][0]][marked[i][1] + 1], TE[i], y);
				else TE[i] = E;
				if (Comp == 0) for (NOff = 0; NOff < Nv; NOff++) Offset[i][NOff] = LTData[0][Iso][NOff][J];
				else for (j = NOff = 0; j < Nv; j++) 
				{
					Offset[i][NOff++] = LTData[0][Iso][j][J-1];
					Offset[i][NOff++] = LTData[0][Iso][j][J+1];
				}
			}
			Spekt->optimizeWeighting(Offset, TE, NOff, NMarked, WF);
			Destroy(Offset, NMarked);
			if ((i = CompBox->currentIndex()) != Comp)
			{
				if (i==0) for (NOff = 0; NOff < Nv; NOff++) WF[NOff] = 0.5 * (WF[2 * NOff] + WF[2 * NOff + 1]);
				else for (NOff = Nv - 1; NOff <= 0; NOff--) WF[2 * NOff] = WF[2 * NOff + 1] = WF[NOff];
				Comp = i;
			}
			Off = new double[NOff];
			v = -1;
		}
		else
		{*/
		TT = UState->getTermTable();
		LT = LState->getTermTable();
		if (LT == 0 && AddByProgFCF->isChecked())
		{
			QMessageBox::information(this, "MolSpektanalysis", "Error: For the lower electronic state no term energies are available!");
			return;
		}
		if (SWF == 0)
		{
			if (TT == 0)
			{
				if (((!AddByBand->isChecked() || AddByProg->isChecked())) || LT == 0)
				{
					QMessageBox::information(this, "MolSpektanalysis", 
						"Error: For at least one state no term energies are available!", QMessageBox::Ok);
					NJ = 0;
					return;
				}
				TData = 0;
				NC = 2;
				Mv = 10000;
			}
			else
			{
				TData = TT->getData();
				NC = TT->getNumComp();
				Mv = TT->getMaxv();
			}
		}
		else Mv = SWF->getNv() - 1;
		FitData *fdat = UState->getFitData();
		QList<LevelComb> LevelList;
		Comp = CompBox->currentIndex();
		if (UpdateAllVs->isChecked())
		{
			if (fdat == 0)
			{
				QMessageBox::information(this, "MolSpektanalysis", 
					"Error: No FitData available!", QMessageBox::Ok);
				NJ = 0;
				return;
			}
			LevelList = fdat->getAvLevelComb(iso, Comp, vmin, vmax);
			if (LogWindow == 0)
			{
				LogWindow = new QTextEdit;
				LogWindow->setReadOnly(true);
				MW->showMDIChild(LogWindow);
			}
		}
		Iso = IsoBox->currentIndex();
		v = L1E->text().toInt();
		vss = L2E->text().toInt();
		if (v > Mv && LevelList.size() == 0)
		{
			QMessageBox::information(this, "MolSpektanalysis", 
				"Error: Vibrational level not available!", QMessageBox::Ok);
			NJ = 0;
			return;
		}
		for (l=0; l < LevelList.size() || (l==0 && !UpdateAllVs->isChecked()); l++)
		{
			if (UpdateAllVs->isChecked())
			{
				Iso = LevelList[l].iso;
				v = LevelList[l].v;
				Comp = 1 - LevelList[l].ef;
			}
			if (v > Mv)
			{
				QMessageBox::information(this, "MolSpektanalysis", 
					"Error: Not all vibrational levels available in term energy table! Please choose the correct one or caculate a new one if the correct one is not available.", 
					QMessageBox::Ok);
				NJ = 0;
				return;
			}
			JStart = UState->getJStart(Iso, Comp);
			JStep = Mol->getJStep(Iso);
			if (SWF == 0)
			{
				MJ = ((TT != 0 ? TT->getMaxJ() < LT->getMaxJ() : false) ? TT->getMaxJ() : LT->getMaxJ());
				c = (1 < NC ? 1 - Comp : 0);
				if (TT != 0) for (J = JStart; (J < MJ ? TData[c][Iso][v][J] != 0 : false); J+= JStep) ;
				else J = 1000;
			}
			else
			{
				ElState *SB = 0;
				int C = -1;
				MJ = (SWF->getNJ() <= LT->getMaxJ() ? SWF->getNJ() - 1 : LT->getMaxJ());
				for (NWFC = c = n = 0; n < SWF->getNChannels(); n++) 
				{
					SWF->getStateComp(n, SB, C);
					if ((TStr[NWFC] = Mol->getTransitionStrength(LState, SB)) > 0.0) WFC[NWFC++] = n;
				}
				I = SWF->getIsoT(Iso);
				if (!SWF->isJA(I, JStart) && SWF->isJA(I, JStart + JStep)) JStart += JStep;
				for (J = JStart; SWF->isJA(I, J) && J < MJ; J += JStep) ;
			}
			if (AddByBand->isChecked())
			{
				LTData = LT->getData();
				for (j = JStart - Comp; (j < MJ ? LTData[0][Iso][v][j] != 0 : false); j += JStep) ;
				if (j + Comp < J) J = j + Comp;
				WF = new double[2];
				WF[0] = WF[1] = 0.5;
				Off = new double[2];
				NOff = Comp + 1;
			}
			else if (AddByProg->isChecked() && TData == 0 && SWF == 0)
			{
				int AM, LMv; 
				Marker *Ma;
				double E = 0.0;
				RSpekt->GetMarker(AM, Ma);
				LTData = LT->getData();
				LMv = LT->getMaxv();
				for (n=m=0; n < AM; n++) if (Ma[n].DisplayData && Ma[n].Iso == RIso && Ma[n].Js == RJs && Ma[n].vs == v)
					if (((Ma[n].Js != Ma[n].Jss && RComp == 1) || (Ma[n].Js == Ma[n].Jss && RComp == 0)) 
						&& (Ma[n].vss <= LMv ? LTData[0][RIso][Ma[n].vss][Ma[n].Jss] != 0.0 : false))
				{
					E += LTData[0][RIso][Ma[n].vss][Ma[n].Jss] + Ma[n].Line[0];
					m++;
				}
				Ra = E/m + Ri;
				Ri = Ra - 2.0 * Ri;
			}
			if (IData != 0) Destroy(IData, NJ);
			IData = Create(NJ = (J - JStart) / JStep, NP); 
			if (SWF != 0)
			{
				JA = new int[NJ];
				UEn = new double[NJ];
				UWF = Create(NJ, NWFC, NumWFPoints);
				vA = new int[NJ];
				CA = new int[NJ];
				for (J = JStart, n=0; n < NJ; n++, J += JStep) JA[n] = J;
				SWF->getWaveFuncs(I, v, WFC, NWFC, NJ, JA, rmin, rmax, NumWFPoints, UWF, UEn, vA, CA);
				delete[] vA;
				vA = 0;
			}
            if (AddByProgFCF->isChecked()) IntensF = new double[NJ];
			for (J = JStart, n=m=0; n < NJ; n++, J += JStep)
			{
				//printf("J'=%d\n", J);
				if (AddByProgFCF->isChecked())
				{
                    if (SWF != 0) Spekt->addByProgFCF(UState, LState, Iso, v, J, Comp, T, Ri, Res, Result, NumWFPoints, UEn[n], UWF[n], TStr, NWFC, IntensF + n);
                    else Spekt->addByProgFCF(UState, LState, Iso, v, J, Comp, T, Ri, Res, Result, NumWFPoints);
				}
				else if (AddByProg->isChecked()) 
				{
					if (TT != 0)
					{
						Ra = TData[c][Iso][v][J] + Rd;
						Ri = TData[c][Iso][v][J] - Rd;
					}
					if (!Spekt->addByProgression(LState, Iso, J, Comp, RSpekt, Ri, Ra, Res, RIso, RJs, RComp, Result)) break;
				}
				else
				{
					if (TT != 0) 
					{
						Ra = TData[c][Iso][v][J] + Rd;
						Ri = TData[c][Iso][v][J] - Rd;
					}
					/*if (Opt) for (i=j=0; i < Nv; i++) for (k = J - Comp; k <= J + Comp; k+=2) Off[j++] = LTData[0][Iso][i][k];
					else*/ for (j=0, k = J - Comp; k <= J + Comp; k+=2) Off[j++] = LTData[0][Iso][vss][k];
					Spekt->add(Off, WF, NOff, Ri, Ra, Res, Result);
				}
				for (i=0; i < NP; i++) IData[n][i] = Result[i][1];
			}
			if (SWF != 0)
			{
				delete[] JA;
				delete[] UWF;
			}
			if (Image == 0) 
			{
				Image = new AddSpectrum(MW, this);
				MW->showMDIChild(Image);
			}
            Image->setData(IData, NJ, JStart, JStep, NP, Spekt, UState, 1 - Comp, Iso, v, TT, Ri, Ra, Res, UEn, vA, CA, IntensF);
			if (RaE->isEnabled()) Image->setRanges(JStart, JStart + (NJ - 1) * JStep, Ri, Ra);
			else Image->setRanges(JStart, JStart + (NJ - 1) * JStep, -Rd, Ri = Rd);
			Image->setUnits("J", "Residuals [cm^{-1}]");
			if (!Image->isVisible()) Image->show();
			Image->activateWindow();
			Image->setFocus();
			if (WF != 0) delete[] WF;
			if (Off != 0) delete[] Off;
			if (UpdateAllVs->isChecked())
			{
				int result = Image->AutoUpdateLevels(true);
				QString LevelText = (LevelList[l].ef == eLevel ? ": e level iso=" :
				                                                  ": f level iso=")
								   + IsoBox->itemText(LevelList[l].iso) 
								   + " v=" + QString::number(LevelList[l].v);
				QString ResultText;
				switch(result)
				{
					case autoUpdateSuccess:
						ResultText = "Success" + LevelText;
						Image->acceptAutoUpdateResult();
						break;
					case autoUpdateInsufficientData:
						ResultText = "Insufficient data" + LevelText;
						break;
					case autoUpdateCompleteFailure:
						ResultText = "Failure" + LevelText;
						break;
					case autoUpdatePartialFailure:
						ResultText = "Partial failure" + LevelText;
						Image->acceptAutoUpdateResult();
						break;
					default:
						ResultText = "Unknown result" + LevelText;
						break;
				}
				LogWindow->append(ResultText);
				printf("%s\n", ResultText.toLatin1().data());
			}
		}
		Destroy(Result, NP);
		Result = 0;
	}
	else
	{
		if (AddByBand->isChecked()) 
			S = Spekt->addByBand(UState, LState, IsoBox->currentIndex(), L1E->text().toInt(),
			    			 L2E->text().toInt(), CompBox->currentIndex(), T, Ri, Res, Result);
		else if (AddByProg->isChecked())
		{
			if (RSpekt == 0)
			{
				QMessageBox::information(this, "MolSpektanalysis", 
						"Error: The selected reference spectrum doesn't exist anymore!", QMessageBox::Ok);
				return;
			}
			S = Spekt->addByProgression(LState, IsoBox->currentIndex(), L1E->text().toInt(),
										CompBox->currentIndex(), RSpekt, Ri, Ra, Res,
										RIsoBox->currentIndex(), RJsE->text().toInt(), 
										RCompBox->currentIndex(), Result);
		}
		else if (AutoAnalyse->isChecked())
			Spekt->autoAddByProgFCF(UState, LState, IsoBox->currentIndex(), CompBox->currentIndex(), T, Ri,
                                    Res, NumWFPoints);
		else
			S = Spekt->addByProgFCF(UState, LState, IsoBox->currentIndex(), L1E->text().toInt(),
                                    L2E->text().toInt(), CompBox->currentIndex(), T, Ri, Res, Result, NumWFPoints);
	}
	if (Result != 0)
	{
		if (S)
		{
			//printf("AddDialog::Apply, n: Result=%d, NP=%d\n", Result, NP);
			Spektrum *nSpektrum = MW->CreateSpectrum();
			if (nSpektrum != 0)
			{
				nSpektrum->setData(Result, NP);
				QString FileName = Spekt->getFileName();
				if ((n = FileName.indexOf(".")) == -1) n = FileName.length();
				nSpektrum->setFileName(FileName.left(n) + "_add" 
										+ FileName.right(FileName.length() - n));
				nSpektrum->setName(FileName.left(n) + "_add");
				nSpektrum->show();
			}
		}
		Destroy(Result, NP);
	}
}
