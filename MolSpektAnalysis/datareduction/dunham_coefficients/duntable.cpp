//
// C++ Implementation: duntable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#include "duntable.h"
#include "molecule.h"
#include "utils.h"
#include "fit.h"
#include "termtable.h"
#include "linetable.h"
#include "elstate.h"
#include "fitdata.h"
#include "isotab.h"
#include "tableline.h"
#include "termenergy.h"

#include <math.h>
#include <limits>

#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QFileDialog>
#include <QStringList>
#include <QGridLayout>

using std::numeric_limits;

DunTable::DunTable(MainWindow *MW, Molecule *Mol) 
	: TableWindow(DunhamTable, MW, Mol)
{
	numCoefficients = vp1 = vp2 = Jp1 = Jp2 = 0;
	State = 0;
	BlockIC = false;
	Tab->setColumnCount(5);
	Tab->setRowCount(MaxDunCoefficients);
	Tab->setHorizontalHeaderLabels(QStringList() << "k" << "l" << "type" << "coefficient" 
			<< "error");
	setWindowTitle("New Dunham coefficient set");
	connect(Tab, SIGNAL(cellChanged(int, int)), this, SLOT(itemChanged(int, int)));
	Saved();
	printf("End DunTable::DunTable\n");
}


DunTable::~DunTable()
{
}

void DunTable::getLinePoints(int &rvp1, int &rJp1, int &rvp2, int &rJp2)
{
	rvp1 = vp1;
	rJp1 = Jp1;
	rvp2 = vp2;
	rJp2 = Jp2;
}

void DunTable::setLinePoints(int svp1, int sJp1, int svp2, int sJp2)
{
	bool changed = (svp1 != vp1 || sJp1 != Jp1 || svp2 != vp2 || sJp2 != Jp2);
	vp1 = svp1;
	Jp1 = sJp1;
	vp2 = svp2;
	Jp2 = sJp2;
	if (changed) Changed();
}

int *DunTable::getmaxv()
{
	int maxJ = getJMax();
	if (maxJ <= 0) maxJ = cMaxJ;
	int *r = new int[maxJ+1], n;
	if ((vp1 == 0 && vp2 == 0) || Jp1 == Jp2)
	{
		int v = getvMax();
		if (v==-1) v = cMaxv;
		for (n=0; n <= maxJ; n++) r[n] = v;
	}
	else
	{
		double g, a, J;
		if (Jp2 < Jp1)
		{
			n = Jp2;
			Jp2 = Jp1;
			Jp1 = n;
			n = vp2;
			vp2 = vp1;
			vp1 = n;
		}
		g = double(vp2 - vp1) / double(Jp2 - Jp1);
		a = double(vp1) - g * double(Jp1);
		for (J=0.0; J <= maxJ; J++) r[int(J)] = int(a + g*J);
	}
	return r;
}

void DunTable::getData(int &kMax, int *&lMax, double **&Par, double **&Korr, double **&SpinRot, double **&adCorr)
{
	int i, j, l;
	bool ldA = false, SpinRA = false, adCorrA = false;
	QString T;
	if (numCoefficients > Tab->rowCount()) numCoefficients = Tab->rowCount();
	for (i=0, kMax = 0; i < numCoefficients; i++)
	{
		if ((j = Tab->item(i, 0)->text().toInt()) > kMax) kMax = j;
		T = Tab->item(i, 2)->text();
		if (T == "ldcorr" || T.toInt() == 6) ldA = true;
		else if (T == "SpinRot") SpinRA = true;
		else if (T == "adCorr") adCorrA = true;
	}
	lMax = new int[kMax+1];	
	for (i=0; i <= kMax; i++) lMax[i] = 0;
	for (i=0; i < numCoefficients; i++) 
		if ((j = Tab->item(i, 1)->text().toInt()) > lMax[l=Tab->item(i,0)->text().toInt()]) 
			lMax[l] = j;
	Par = new double*[kMax+1];
	Korr = (ldA ? new double*[kMax+1] : 0);
	SpinRot = (SpinRA ? new double*[kMax+1] : 0);
	adCorr = (adCorrA ? new double*[kMax+1] : 0);
	for (i=0; i <= kMax; i++)
	{
		for (j=0, Par[i] = new double[lMax[i] + 1]; j <= lMax[i]; j++) Par[i][j] = 0.0; 
		if (ldA) for (j=0, Korr[i] = new double[lMax[i] + 1]; j <= lMax[i]; j++) Korr[i][j] = 0.0;
		if (SpinRA) for (j=0, SpinRot[i] = new double[lMax[i] + 1]; j <= lMax[i]; j++) SpinRot[i][j] = 0.0;
		if (adCorrA) for (j=0, adCorr[i] = new double[lMax[i] + 1]; j <= lMax[i]; j++) adCorr[i][j] = 0.0;
	}
	for (i=0; i < numCoefficients; i++) 
	{
		T = Tab->item(i, 2)->text();
		if ((j = T.toInt()) == 1 || T == "coeff") 
			Par[Tab->item(i, 0)->text().toInt()][Tab->item(i, 1)->text().toInt()] 
												= Tab->item(i, 3)->text().toDouble();
		else if (j == 6 || T == "ldcorr") 
			Korr[Tab->item(i, 0)->text().toInt()][Tab->item(i, 1)->text().toInt()]
					= Tab->item(i, 3)->text().toDouble();
		else if (T == "SpinRot")
			SpinRot[Tab->item(i, 0)->text().toInt()][Tab->item(i, 1)->text().toInt()] = Tab->item(i, 3)->text().toDouble();
		else if (T == "adCorr")
			adCorr[Tab->item(i, 0)->text().toInt()][Tab->item(i, 1)->text().toInt()] = Tab->item(i, 3)->text().toDouble();
	}
}

void DunTable::getBandConstants(double **&C, double **&err, int &NC, int &MLD, int &MSR, int &MAD, int &nv)
{
	int k, l, Mv = getvMax(), kmax, *lmax, MC = -1, RC = Tab->rowCount(), v;
	double **DK, **LDK, **SR, **ADK, **DKerr, **LDKerr = 0, **SRerr = 0, **ADKerr = 0, vF, vP;
	QString T;
	if (Mv <= 0) Mv = cMaxv;
	getData(kmax, lmax, DK, LDK, SR, ADK);
	DKerr = new double*[kmax + 1];
	for (k=0; k <= kmax; k++)
	{
		if (lmax[k] > MC) MC = lmax[k];
		DKerr[k] = new double[lmax[k] + 1];
	}
	NC = (++MC);
	if (LDK != 0)
	{
		LDKerr = new double*[kmax + 1];
		for (k = MLD = 0; k <= kmax; k++) 
		{
			for (l=0; l <= lmax[k]; l++) if (LDK[k][l] != 0.0 && l > MLD) MLD = l;
			LDKerr[k] = new double[lmax[k] + 1];
		}
		NC += (++MLD);
	}
	else MLD = 0;
	if (SR != 0)
	{
		SRerr = new double*[kmax + 1];
		for (k = MSR = 0; k <= kmax; k++) 
		{
			for (l=0; l <= lmax[k]; l++) if (SR[k][l] != 0.0 && l > MSR) MSR = l;
			SRerr[k] = new double[lmax[k] + 1];
		}
		NC += (++MSR);
	}
	else MSR = 0;
	if (ADK != 0)
	{
		ADKerr = new double*[kmax + 1];
		for (k = MAD = 0; k <= kmax; k++) 
		{
			for (l=0; l <= lmax[k]; l++) if (ADK[k][l] != 0.0 && l > MAD) MAD = l;
			ADKerr[k] = new double[lmax[k] + 1];
		}
		NC += (++MAD);
	}
	else MAD = 0;
	C = Create(nv = Mv + 1, NC);
	err = Create(nv, NC);
	for (k=0; k < RC; k++)
	{
		T = Tab->item(k, 2)->text();
		if ((l = T.toInt()) == 1 || T == "coeff") 
			DKerr[Tab->item(k, 0)->text().toInt()][Tab->item(k, 1)->text().toInt()] = Tab->item(k, 4)->text().toDouble();
		else if (l == 6 || T == "ldcorr") 
			LDKerr[Tab->item(k, 0)->text().toInt()][Tab->item(k, 1)->text().toInt()] = Tab->item(k, 4)->text().toDouble();
		else if (T == "SpinRot")
			SRerr[Tab->item(k, 0)->text().toInt()][Tab->item(k, 1)->text().toInt()] = Tab->item(k, 4)->text().toDouble();
		else if (T == "adCorr")
			ADKerr[Tab->item(k, 0)->text().toInt()][Tab->item(k, 1)->text().toInt()] = Tab->item(k, 4)->text().toDouble();
	}
	for (k=0; k <= kmax; k++) for (l=0; l <= lmax[k]; l++)
	{
		DKerr[k][l] *= DKerr[k][l];
		if (LDK != 0) LDKerr[k][l] *= LDKerr[k][l];
		if (SR != 0) SRerr[k][l] *= SRerr[k][l];
		if (ADK != 0) ADKerr[k][l] *= ADKerr[k][l];
	}
	for (v=0; v < nv; v++) for (l=0; l < NC; l++)
	{
		C[v][l] = 0.0;
		err[v][l] = 0.0;
	}
	for (v=0, vF = 0.5; v < nv; v++, vF += 1.0) for (k=0, vP = 1.0; k <= kmax; k++, vP *= vF)
	{
		for (l=0; l < MC && l <= lmax[k]; l++)
		{
			C[v][l] += vP * DK[k][l];
			err[v][l] += vP * vP * DKerr[k][l];
		}
		for (l=0; l < MLD && l <= lmax[k]; l++)
		{
			C[v][MC + l] += vP * LDK[k][l];
			err[v][MC + l] += vP * vP * LDKerr[k][l];
		}
		for (l=0; l < MSR && l <= lmax[k]; l++)
		{
			C[v][MC + MLD + l] += vP * SR[k][l];
			err[v][MC + MLD + l] += vP * vP * SRerr[k][l];
		}
		for (l=0; l < MAD && l <= lmax[k]; l++)
		{
			C[v][MC + MLD + MSR + l] += vP * ADK[k][l];
			err[v][MC + MLD + MSR + l] += vP * vP * ADKerr[k][l];
		}
	}
	if (MC > 2) for (v=0; v < nv; v++) C[v][2] *= -1;
	for (v=0; v < nv; v++) for (l=0; l < NC; l++) err[v][l] = sqrt(err[v][l]);
	Destroy(DK, kmax + 1);
	Destroy(DKerr, kmax + 1);
	if (LDK != 0)
	{
		Destroy(LDK, kmax + 1);
		Destroy(LDKerr, kmax + 1);
	}
	if (SR != 0)
	{
		Destroy(SR, kmax + 1);
		Destroy(SRerr, kmax + 1);
	}
	if (ADK != 0)
	{
		Destroy(ADK, kmax + 1);
		Destroy(ADKerr, kmax + 1);
	}
}

double DunTable::getwe()
{
	int i;
	for (i=0; i < numCoefficients; i++) 
		if (Tab->item(i, 0)->text().toInt() == 1 && Tab->item(i, 1)->text().toInt() == 0
			&& (Tab->item(i, 2)->text().toInt() == 1 || Tab->item(i, 2)->text() == "coeff"))
				return Tab->item(i, 3)->text().toDouble();
	return 0.0;
}

void DunTable::setData(int N, int *t, int *k, int *l, double *C)
{
	int n;
	BlockIC = true;
	Tab->setRowCount(numCoefficients = N);
	for (n=0; n<N; n++)
	{
		Tab->setItem(n, 0, new QTableWidgetItem(QString::number(k[n])));
		Tab->setItem(n, 1, new QTableWidgetItem(QString::number(l[n])));
		switch (t[n])
		{
			case 1:
				Tab->setItem(n, 2, new QTableWidgetItem("coeff"));
				break;
			case 6:
				Tab->setItem(n, 2, new QTableWidgetItem("ldcorr"));
				break;
			default:
				Tab->setItem(n, 2, new QTableWidgetItem(QString::number(t[n])));
				break;
		}
		Tab->setItem(n, 3, new QTableWidgetItem(QString::number(C[n], 'g', 16)));
		Tab->setItem(n, 4, new QTableWidgetItem(""));
	}
	BlockIC = false;
}

void DunTable::addTableLine()
{
	Tab->setRowCount(Tab->rowCount() + 1);
}

void DunTable::calcTermEnergies(TermTable *&TT, bool show, int vMax, int JMax)
{
	if (numCoefficients == 0)
	{
		printf("DunTable::calcTermEnergies: Error: no coefficients available!\n");
		return;
	}
	if (State == 0)
	{
		printf("DunTable::calcTermEnergies: Error: no ElState available!\n");
		return;
	}
	if (molecule == 0)
	{
		printf("DunTable::calcTermEnergies: Error: no molecule available!\n");
		return;
	}
	if (MW == 0)	
	{
		printf("DunTable::calcTermEnergies: Error: no MainWindow available!\n");
		return;
	}
	
	int i, j, k, l, J, v, *lMax, kMax, NumIso, NumComp, Nef, NFC, ef, FC, c, lM;
	if (vMax < 0) 
	{
		vMax = getvMax();
		if (vMax <= 0) vMax = cMaxv;
	}
	if (JMax < 0) 
	{
		JMax = getJMax();
		if (JMax <= 0) JMax = cMaxJ;
	}
	int *mv = getmaxv();
	double T;
	IsoTab *Iso = molecule->getIso();
	//printf("l=%d\n", l);
	k = State->getOmega();
	int OmegaSQ = k*k;
	NumIso = Iso->numIso;
	double **DKoeff, **DKorr, **SpinRot, ****SpinRF, **adCorr;
	getData(kMax, lMax, DKoeff, DKorr, SpinRot, adCorr);
	NumComp = (Nef = (DKorr != 0 ? 2 : 1)) * (NFC = (SpinRot != 0 ? 2 : 1));
	//for (k=0; k <= kMax; k++) for (l=0; l <= lMax[k]; l++) 
		//	printf("k=%d, l=%d, DKoeff=%f, DKorr=%f\n", k, l, DKoeff[k][l], DKorr[k][l]);
	double ****TermData = Create(NumComp, NumIso, vMax + 1, JMax + 1); 
	double **vF = Create(NumIso, vMax + 1);
	double **JF = Create(NumIso, JMax + 1);
	if (SpinRot != 0)
	{
		for (k=1, lM = lMax[0]; k <= kMax; k++) if (lMax[k] > lM) lM = lMax[k];
		SpinRF = cSpinRotJF(Iso, JMax, lM);
	}
	else SpinRF = 0;
	for (i=0; i<NumIso; i++) 
	{
		//printf("i=%d\n", i);
		vF[i][0] = Iso->rootRRM[i] * 0.5;
		for (J=0; J<=JMax; J++) for (j=0; j < NumComp; j++) TermData[j][i][0][J] = 0.0;
		for (v=1; v<=vMax; v++) 
		{
			vF[i][v] = vF[i][v-1] + Iso->rootRRM[i];
			for (J=0; J<=JMax; J++) for (j=0; j < NumComp; j++) TermData[j][i][v][J] = 0.0;
		}
		JF[i][0] = -1.0 * OmegaSQ * Iso->relRedMass[i];
		for (J=1; J<=JMax; J++) JF[i][J] = JF[i][J-1] + 2.0 * double(J) * Iso->relRedMass[i];
	}
	for (ef = 0; ef < Nef; ef++) for (FC = 0; FC < NFC; FC++) 
		for (k = kMax, c = Nef * FC + ef; k >= 0; k--) for (v=0; v<=vMax; v++) for (i=0; i<NumIso; i++) 	
	{
		for (J=0; J<=JMax; J++)
		{
			TermData[c][i][v][J] *= vF[i][v];
			T = 0.0;
			for (l=lMax[k]; l>=0; l--)
			{
				T *= JF[i][J];
				T += (DKorr != 0 && ef == 0 ? DKoeff[k][l] + DKorr[k][l] : DKoeff[k][l]);
				if (adCorr != 0 ? adCorr[k][l] != 0.0 : false) T += (1.0 - Iso->relRedMass[i]) * adCorr[k][l];
			}
			TermData[c][i][v][J] += T;
			if (SpinRot != 0) for (l=0; l <= lMax[k]; l++) TermData[c][i][v][J] += SpinRF[FC][i][J][l] * SpinRot[k][l];
			//if (i==0 && v==0 && J==0) printf("k=%d, T=%f\n", k, T);
		}
	}
	if (mv != 0) for (J=0; J<=JMax; J++) for (v = (mv[J] > -1 ? mv[J] + 1 : 0); v <= vMax; v++) 
				for (i=0; i<NumIso; i++) for (j=0; j < NumComp; j++) TermData[j][i][v][J] = 0.0;
	for (i=0; i<= kMax; i++)
	{
		delete[] DKoeff[i];
		if (DKorr != 0) delete[] DKorr[i];
		if (SpinRot != 0) delete[] SpinRot[i];
		if (adCorr != 0) delete[] adCorr[i];
	}
	delete[] DKoeff;
	if (DKorr != 0) delete[] DKorr;
	if (SpinRot != 0) delete[] SpinRot;
	if (adCorr != 0) delete[] adCorr;
	delete[] lMax;
	//for (J=0; J<=JMax; J++) printf("TermData[0][0][0][%d]=%f\n", J, TermData[0][0][0][J]);
	Destroy(vF, Iso->numIso);
	Destroy(JF, Iso->numIso);
	if (SpinRF != 0) Destroy(SpinRF, 2, Iso->numIso, JMax + 1);
	delete Iso;
	if (TT == 0) 
	{
		TT = MW->CreateTermTable();
		j = molecule->getNumStates();
		for (i=0; i < j; i++) if (molecule->getDunK(i) == this) molecule->addTerm(i, TT);
		TT->setName("Term" + getName());
		TT->setSource("Calculated from " + getSource());
	}
	TermTable *XTT = molecule->getTerm(0);
	int *CompZ = 0;
	if (XTT != 0 ? XTT->getNumComp() == NumComp : false)
	{
		CompZ = new int[NumComp];
		int *XCompZ = XTT->getCompZ();
		for (i=0; i < NumComp; i++) CompZ[i] = XCompZ[i];
	}
	TT->setData(TermData, NumComp, NumIso, vMax, JMax, CompZ);
	if (show) TT->show();
}

bool DunTable::isLambdaDoublingAv()
{
	int n;
	if (numCoefficients > Tab->rowCount()) numCoefficients = Tab->rowCount();
	for (n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "ldcorr" || Tab->item(n, 2)->text().toInt() == 6)
		return true;
	return false;
}

bool DunTable::isSpinRAv()
{
	int n;
	if (numCoefficients > Tab->rowCount()) numCoefficients = Tab->rowCount();
	for (n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "SpinRot") 
		return true;
	return false;
}

void DunTable::removeFinestructure(TableLine *&Data, int &N)
{
	FitData *FD = (State != 0 ? State->getFitData() : 0);
	if (FD == 0) 
	{
		N=0;
		Data = 0;
		return;
	}
	int n, mv = 0, mJ = 0, mvF = 0, mJF = 0, i, nSR = 0;
	IsoTab *Iso = molecule->getIso();
	FD->getData(Data, N);
	for (n=0; n<N; n++) 
	{
		if (Data[n].vss > mv) mv = Data[n].vss;
		if (Data[n].Jss > mJ) mJ = Data[n].Jss;
	}
	for (n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "SpinRot")
	{
		nSR++;
		if ((i = Tab->item(n, 0)->text().toInt()) > mvF) mvF = i;
		if ((i = Tab->item(n, 1)->text().toInt()) > mJF) mJF = i;
	}
	double ***vF = cvF(Iso->numIso, Iso, mv, mvF), ****JF = cSpinRotJF(Iso, mJ, mJF);
	double SR[nSR];
	int SRvF[nSR], SRJF[nSR];
	for (n=i=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "SpinRot")
	{
		SRvF[i] = Tab->item(n, 0)->text().toInt();
		SRJF[i] = Tab->item(n, 1)->text().toInt();
		SR[i++] = Tab->item(n, 3)->text().toDouble();
	}
	for (n=0; n<N; n++) 
	{	
		if (Data[n].FC >= 0 && Data[n].FC <= 1) for (i=0; i < nSR; i++)
			Data[n].WN += (Data[n].isTE ? -SR[i] : SR[i])
							* vF[Data[n].Iso][Data[n].vss][SRvF[i]]
							* JF[Data[n].FC][Data[n].Iso][Data[n].Jss][SRJF[i]];
	}
	Destroy(vF, Iso->numIso, mv + 1);
	Destroy(JF, 2, Iso->numIso, mJ + 1);
	delete Iso;
}

void DunTable::removeLambdaDoubling(TableLine *&Data, int &N)
{
	FitData *FD = (State!= 0 ? State->getFitData() : 0);
	if (FD == 0) return;
	int n, mv = 0, mJ = 0, mvF = 0, mJF = 0, i, j, nLD = 0;
	IsoTab *Iso = molecule->getIso();
	FD->getData(Data, N);
	for (n=0; n<N; n++)
	{
		if (Data[n].vss > mv) mv = Data[n].vss;
		if (Data[n].Jss > mJ) mJ = Data[n].Jss;
	}
	for (n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "ldcorr" || Tab->item(n, 2)->text().toInt() == 6)
	{
		nLD++;
		if ((i = Tab->item(n, 0)->text().toInt()) > mvF) mvF = i;
		if ((i = Tab->item(n, 1)->text().toInt()) > mJF) mJF = i;
	}
	double ***vF = cvF(Iso->numIso, Iso, mv, mvF), ***JF = cJF(Iso->numIso, Iso, mJ, mJF, State->getOmega());
	double LD[nLD];
	int LDvF[nLD], LDJF[nLD];
	for (n=i=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "ldcorr" || Tab->item(n, 2)->text().toInt() == 6)
	{
		LDvF[i] = Tab->item(n, 0)->text().toInt();
		LDJF[i] = Tab->item(n, 1)->text().toInt();
		LD[i++] = Tab->item(n, 3)->text().toDouble();
	}
	for (n=0; n<N; n++) if (Data[n].Js != Data[n].Jss) for (j=0; j < nLD; j++) 
		Data[n].WN -= LD[j] * vF[Data[n].Iso][Data[n].vss][LDvF[j]] * JF[Data[n].Iso][Data[n].Jss][LDJF[j]];
	Destroy(vF, Iso->numIso, mv + 1);
	Destroy(JF, Iso->numIso, mJ + 1);
	delete Iso;
}

void DunTable::applyAdiabaticCorrection(TableLine*& Data, int& N)
{
	FitData *FD = (State != 0 ? State->getFitData() : 0);
	if (FD == 0) return;
	int n, mv = 0, mJ = 0, mvF = 0, mJF = 0, i, j, nAC = 0;
	IsoTab *Iso = molecule->getIso();
	FD->getData(Data, N);
	for (n=0; n<N; n++)
	{
		if (Data[n].vss > mv) mv = Data[n].vss;
		if (Data[n].Jss > mJ) mJ = Data[n].Jss;
	}
	for (n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "adCorr")
	{
		nAC++;
		if ((i = Tab->item(n, 0)->text().toInt()) > mvF) mvF = i;
		if ((i = Tab->item(n, 1)->text().toInt()) > mJF) mJF = i;
	}
	double ***vF = cvF(Iso->numIso, Iso, mv, mvF), ***JF = cJF(Iso->numIso, Iso, mJ, mJF, State->getOmega()), AC[nAC];
	int ACvF[nAC], ACJF[nAC];
	for (i=n=0; n < numCoefficients; n++) if (Tab->item(n, 2)->text() == "adCorr")
	{
		ACvF[i] = Tab->item(n, 0)->text().toInt();
		ACJF[i] = Tab->item(n, 1)->text().toInt();
		AC[i++] = Tab->item(n, 3)->text().toDouble();
	}
	for (n=0; n<N; n++) for (j=0; j < nAC; j++) 
		Data[n].WN += AC[j] * (Iso->relRedMass[Data[n].Iso] - 1.0) * vF[Data[n].Iso][Data[n].vss][ACvF[j]]
							* JF[Data[n].Iso][Data[n].Jss][ACJF[j]];
	Destroy(vF, Iso->numIso, mv + 1);
	Destroy(JF, Iso->numIso, mJ + 1);
	delete Iso;
}

bool DunTable::readData(QString FileName)
{
	//printf("Beginn DunTable::readData\n");
	int n, m, i, j, k, l, r = 0;
	bool nT = true, R = false;
	QString Buffer;
	QStringList L;
	QFile Datei(FileName);
	if (!read(&Datei)) return false;
	QTextStream S(&Datei);
	Buffer = S.readLine();
	if (Buffer.left(7) != "Source:" && Buffer.left(5) != "Name:" && !Buffer.isEmpty())
	{
		BlockIC = true;
		setSource(Buffer);
		setImported();
		setName("ImportedDunhamCoefficients");
		int sn=-1, asn, N = S.readLine().mid(5, 5).toInt();
		QString T;
		for (n=0; n<N; n++)
		{
			Buffer = S.readLine();
			Tab->setItem(r, 0, new QTableWidgetItem(Buffer.left(5)));
			Tab->setItem(r, 1, new QTableWidgetItem(Buffer.mid(5, 5)));
			if ((T = Buffer.mid(10, 5)) == "    1") Tab->setItem(r, 2, new QTableWidgetItem("coeff"));
			if ((asn = Buffer.mid(15, 5).toInt()) != sn)
			{
				if (sn == -1) sn = asn;
				else break;
			}
			Tab->setItem(r, 3, new QTableWidgetItem(Buffer.mid(20, 22)));
			Tab->setItem(r++, 4, new QTableWidgetItem(""));
		}
		numCoefficients = r;
		if (n<N && MW != 0)
		{
			int NC = N-r;
			int k[NC], l[NC], t[NC];
			double K[NC];
			for (r=0; r < NC; r++)
			{
				if (r>0) Buffer = S.readLine();
				k[r] = Buffer.left(5).toInt();
				l[r] = Buffer.mid(5, 5).toInt();
				t[r] = Buffer.mid(10, 5).toInt();
				K[r] = Buffer.mid(20, 22).toDouble();
			}
			DunTable *Tab2 = MW->CreateDunTable();
			Tab2->setSource(getSource());
			Tab2->setName("ImportedDunhamCoefficients");
			Tab2->setData(NC, t, k, l, K);
			Tab2->show();
		}
		r = numCoefficients;
		BlockIC = false;
	}
	else if (Buffer.left(3) == "DTZ")
	{
		//printf("Begin readData sec DTZ\n");
		Buffer = S.readLine();
		int lm = Buffer.right(Buffer.length() - Buffer.indexOf("=") - 1).toInt();
		int km[lm + 1];
		Buffer = S.readLine();
		L = Buffer.right(Buffer.length() - Buffer.indexOf("=") - 1).split("	",
                    Qt::SkipEmptyParts);
		if (L.size() <= lm) 
		{
			printf("DunTable::readData: error L.size()=%d <= lm=%d\n", L.size(), lm);
			QMessageBox::warning(this, tr("QT4MolSpektAn"), 
							 "Error while opening file " + FileName + '!');
			return false;
		}
		for (n=0, m=0; n <= lm; n++) if ((km[n] = L[n].toInt()) > m) m = km[n];
		//printf("Mitte\n");
		BlockIC = true;
		for (n=0; n <= m; n++)
		{
            L = S.readLine().split("	", Qt::SkipEmptyParts);
			l = L.size();
			for (i=0, k=-1; i < l; i++)
			{
				for (j=0; L[i][j] == ' '; j++) ;
				if (j == L[i].length()) continue; 
				if ((j = L[i].indexOf(",")) != -1) L[i][j] = '.';
				while (k < lm ? n > km[++k] : false) ;
				Tab->setItem(r, 0, new QTableWidgetItem(QString::number(n)));
				Tab->setItem(r, 1, new QTableWidgetItem(QString::number(k)));
				Tab->setItem(r, 2, new QTableWidgetItem("coeff"));
				Tab->setItem(r++, 3, new QTableWidgetItem(L[i]));
			}
		}
		BlockIC = false;
		//printf("Nach Schleife\n");
        L = S.readLine().split("	", Qt::SkipEmptyParts);
		setImported();
		if (L.size() < 2) 
		{
			printf("DunTable::readData: error L.size()=%d < 2\n", L.size());
			QMessageBox::warning(this, tr("QT4MolSpektAn"), 
							 "Error while opening file " + FileName + '!');
			return false;
		} 
		//printf("vMax=%s\n", vMax->text().ascii());
		vMax->setText(L[0].right(L[0].length() - L[0].indexOf("=") - 1));
		//printf("Zeile vorher\n");
		JMax->setText(L[1].right(L[1].length() - L[1].indexOf("=") - 1));
		//printf("End readData sec DTZ\n");
	}
	else
	{		
		printf("NF1\n");
		BlockIC = true;
		while (!S.atEnd())
		{
			if (R) Buffer = S.readLine();
			else R = true;
			if (nT)
			{
				if (Buffer.indexOf("Source", 0, Qt::CaseInsensitive) != -1)
				{
					printf("Source\n");
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					setSource(Buffer.right(Buffer.length() - n));
				}
				else if (Buffer.indexOf("Name", 0, Qt::CaseInsensitive) != -1)
				{
					printf("Name\n");
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					setName(Buffer.right(Buffer.length() - n));
				}
				else if (Buffer.indexOf("Max v", 0, Qt::CaseInsensitive) != -1)
				{
					printf("Maxv\n");
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					vMax->setText(Buffer.right(Buffer.length() - n));
				}
				else if (Buffer.indexOf("Max J", 0, Qt::CaseInsensitive) != -1)
				{
					printf("MaxJ\n");
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					JMax->setText(Buffer.right(Buffer.length() - n));
				}
				else if (Buffer.indexOf("Error", 0, Qt::CaseInsensitive) != -1)
				{
					printf("Error\n");
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					error->setText(Buffer.right(Buffer.length() - n));
				}
				else if (Buffer.indexOf("Border line points", 0, Qt::CaseInsensitive) != -1)
				{
					n = Buffer.indexOf("vp1", 0, Qt::CaseInsensitive);
					n = Buffer.indexOf("=", n+1);
					m = Buffer.indexOf(",", n+1);
					if (m == -1) m = Buffer.length();
					vp1 = Buffer.mid(n+1, m-n-1).toInt();
					n = Buffer.indexOf("Jp1", 0, Qt::CaseInsensitive);
					n = Buffer.indexOf("=", n+1);
					m = Buffer.indexOf(",", n+1);
					if (m == -1) m = Buffer.length();
					Jp1 = Buffer.mid(n+1, m-n-1).toInt();
					n = Buffer.indexOf("vp2", 0, Qt::CaseInsensitive);
					n = Buffer.indexOf("=", n+1);
					m = Buffer.indexOf(",", n+1);
					if (m == -1) m = Buffer.length();
					vp2 = Buffer.mid(n+1, m-n-1).toInt();
					n = Buffer.indexOf("Jp2", 0, Qt::CaseInsensitive);
					n = Buffer.indexOf("=", n+1);
					m = Buffer.indexOf(",", n+1);
					if (m == -1) m = Buffer.length();
					Jp2 = Buffer.mid(n+1, m-n-1).toInt();
				}
				else if (Buffer.indexOf("Table", 0, Qt::CaseInsensitive) != -1) nT = false;
			}
			else if (Buffer.indexOf(':') == -1 && r < MaxDunCoefficients && !Buffer.isEmpty())
			{
				L = Buffer.split("|");
				if (L.size() < 4) continue;
				for (i=0; i<4; i++) Tab->setItem(r, i, new QTableWidgetItem(L[i].remove(" ")));
				if (L.size() > 4) Tab->setItem(r, i, new QTableWidgetItem(L[i].remove(" ")));
				r++;
				printf("r=%d\n", r);
			}
		}
		BlockIC = false;
		Saved();
	}
	Tab->setRowCount(r);
	numCoefficients = r;
	return true;
}

void DunTable::setElState(ElState *elState)
{
	State = elState;
}

bool DunTable::writeData(QString NFilename)
{
	printf("DunTable::writeData\n");
	int i, j;
	QString N = getName();
	QFile Datei(NFilename);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
	S << "Source: " << getSource() << "\n";
	S << "Name: " << getName() << "\n";
	S << "Max v: " << vMax->text() << "\n";
	S << "Max J: " << JMax->text() << "\n";
	S << "Error: " << error->text() << "\n";
	S << "Border line points: vp1 = " << QString::number(vp1) << ", Jp1 = " << QString::number(Jp1) << ", vp2 = " << QString::number(vp2) 
			<< ", Jp2 = " << QString::number(Jp2) << "\n\n";
	S << "Table of coefficients: " << "\n";
	S << "k: | l: | type: | coefficient: | error:\n";
	for (i=0; i < numCoefficients; i++)
	{
		for (j=0; j < 3; j++) S << Tab->item(i, j)->text() << " | ";
		S << Tab->item(i, 3)->text();
		if (Tab->columnCount() > 4) S << " | " << Tab->item(i, 4)->text();
		S << "\n";
		//printf("i=%d\n", i);
	}
	return true;
}

void DunTable::exportTF()
{
	QString FN = QFileDialog::getSaveFileName(this, tr("Export data"), "", 
											  tr("All file names (*.*)"));
	QFile Datei(FN);
	if (!Datei.open(QIODevice::WriteOnly | QIODevice::Truncate))
	{
		QMessageBox::warning(this, tr("QT4MolSpektAn"), 
						 tr("Error while opening/creating file!"));
		return;
	}
	int r, R = Tab->rowCount();
	QTextStream S(&Datei);
	S << getName() << ", " << getSource() << "\n";
	S << "    0" << ("     " + QString::number(R)).right(5) << "    1    1    0    1.0\n";
	for (r=0; r<R; r++) 
		S << ("     " + Tab->item(r, 0)->text()).right(5) 
				<< ("     " + Tab->item(r, 1)->text()).right(5) << "    1    1"
				<< ("                       " + Tab->item(r, 3)->text()).right(22) << "    0    2\n";
}

QString DunTable::getTexTable()
{
	QStringList L;
	int i, n, m, p=0, ml, Ml, l, Maxl=0, k, Maxk=0, Mk, nC, ND, NP, mv, mJ;
	QString Buffer, NB;
	bool FS;
	FitData *FD = (State != 0 ? State->getFitData() : 0);
	if (FD == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"To use this function, the table of Dunham coefficients has to be assigned to the electronic state where the fitdata set is selected, the coefficients are fitted to!");
		return "";
		
	}
	//Transition *ES[MaxStates], *GS[MaxStates];
	int **DunList =  getDunList(nC, n, FS, Maxk, Maxl), *RowN = new int[FD->getNumLines()];
	double **CS = Create(Maxk + 1, Maxl + 1), **Err=0, mp, B;
	int **nD = CreateInt(Maxk + 1, Maxl + 1), lang, nCol;
	QDialog *D = new QDialog(this);
	QGridLayout *La = new QGridLayout(D);
	D->setWindowTitle("Select language and number of columns");
	La->addWidget(new QLabel("Number of colums:", D), 0, 0);
	QComboBox *ColB = new QComboBox(D);
	ColB->setEditable(false);
	ColB->addItems(QStringList() << "2" << "3" << "4" << "5" << "6" << "7" << "8");
	La->addWidget(ColB, 0, 1);
	La->addWidget(new QLabel("Language:", D), 1, 0);
	QComboBox *LangB = new QComboBox(D);
	LangB->setEditable(false);
	LangB->addItems(QStringList() << "english" << "german");
	La->addWidget(LangB, 1, 1);
	La->setRowMinimumHeight(2, 20);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	La->addWidget(OK, 3, 0);
	La->addWidget(Cancel, 3, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Rejected)
	{
		delete D;
		return "";
	}
	lang = LangB->currentIndex();
	nCol = ColB->currentIndex() + 2;
	delete D;
	for (l=0; l <= Maxl; l++) for (k=0; k <= Maxk; k++) 
	{
		CS[k][l] = 0.0;
		nD[k][l] = 0;
	}
	mJ = FD->getMaxJ();
	mv = FD->getMaxv();
	if ((n = getJMax()) < mJ && n != -1) mJ = n;
	else setJMax(mJ);
	int *Mv = getmaxv();
	if (Mv[0] < mv) mv = Mv[0];
	else if (Mv[0] > mv) setvMax(mv);
	for (n=0; n <= mJ; n++) Mv[n] = mv;
	ND = FD->getNumLines(Mv, mJ);
	NP = FD->getNumProgressions(Mv, mJ);
	if (ND > 0)
	{
		IsoTab *Iso = molecule->getIso();
		int xx[ND][5], numD;
		double yy[ND], ssig[ND], SS[NP], Res[nC], err[nC], Tic, dev[ND], V[nC], OV;
		double ***JF = cJF(Iso->numIso, Iso, mJ, Maxl, State->getLambda());
		double ***vF = cvF(Iso->numIso, Iso, mv, Maxk);
		//calcXY(GS, gC, ES, eC, NP, xx, yy, ssig, SS);
		calcXY(FD, NP, ND, xx, yy, ssig, RowN, SS);
		//double **EQS = cEQS(ND, nC, DunList, vF, JF, xx, yy);
		double FQS, aFQS = fitDun(ND, nC, DunList, vF, JF, Res, err, xx, yy, ssig, SS, dev, Iso);
		//printf("Original FQS=%f\n", aFQS);
		//printf("nC=%d\n", nC);
		for (i=0; i < nC; i++) V[i] = 0.0;
		for (i=0, FQS = aFQS; i < numCoefficients; i++)
		{
			for (n=0, mp = 0.0; n < nC - i; n++) if ((B = fabs(err[n] / Res[n])) > mp)
			{
				mp = B;
				p=n;
			}
			for (n=0, m=-1; m<p; n++) if (V[n] == 0.0) m++;
			n--;
			numD = int(floor(log10(err[p])));
			//printf("i=%d\n", i);
			if (i==0) numD--;
			Tic = pow(10.0, numD);
			numD = int(floor(log10(fabs(Res[p])))) - numD;
			for (FQS *= 2.0, OV = Res[p]; FQS > 1.01 * aFQS; Tic /= 10.0, numD++)
			{
				B = Tic * floor(OV / Tic);
				//printf("Res[%d]=%g, err[%d]=%g, Tic[%d]=%g", i, Res[i], i, err[i], i, Tic[i]);
				V[n] = (OV - B > 0.5 * Tic ? B + Tic : B);
				//printf(" Res[%d,r]=%g\n", i, Res[i]);
				FQS = fitDun(ND, nC, DunList, vF, JF, Res, err, xx, yy, ssig, SS, dev, Iso, V);
				printf("V[%d]=%e, numD=%d, FQS=%f\n", n, V[n], numD, FQS);
			}
			if (DunList[n][2] == 0)
			{
				CS[DunList[n][0]][DunList[n][1]] = V[n];
				nD[DunList[n][0]][DunList[n][1]] = numD;
			}
			Tab->item(n, 3)->setText(QString::number(V[n], 'g', numD + 1));
			FQS = aFQS;
		}
		Destroy(vF, Iso->numIso, mv + 1);
		Destroy(JF, Iso->numIso, mJ + 1);
		Changed();
		QMessageBox::information(this, "MolSpektAnalysis", 
							"The final sigma is " + QString::number(sqrt(FQS / double(ND - nC)), 'f', 2));
	}
	else for (i=0; i < numCoefficients; i++) CS[DunList[i][0]][DunList[i][1]]
					= Tab->item(i, 3)->text().toDouble();
	/*if (CS[0][0] > 10.0) 
	{
		T = CS[0][0];
		CS[0][0] = 0.0;
	}*/
	L << "\\begin{table}" << "\\caption{???}";
	L << "\\label{" + (State != 0 ? State->getName() : QString("???"))  + "Dunham}";
	L << "\\centering";
	Buffer = "\\begin{tabular}{r";
	for (n=0; n < nCol; n++) Buffer += 'r';
	L << Buffer + '}';
	for (ml = Ml = 0; Ml <= Maxl; ml = Ml)
	{
		Ml = ml + nCol;
		for (Mk = 0, l = ml; l < Ml && l <= Maxl; l++)
		{
			for (k = Maxk; (k >= 0 ? CS[k][l] == 0.0 : false); k--) ;
			if (k > Mk) Mk = k;
			//printf("l=%d, k=%d, Mk=%d\n", l, k, Mk);
		}
		L << "\\hline\\noalign{\\smallskip}";
		Buffer = "\\(l\\downarrow{} k\\rightarrow{}\\)";
		for (l = ml; l < Ml; l++) Buffer += " & " + (l <= Maxl ? QString::number(l) : QString(" "));
		L << Buffer + "\\\\ \\hline\\noalign{\\smallskip}";
		for (k=0; k <= Mk; k++) 
		{
			Buffer = QString::number(k);
			for (l = ml; l < Ml; l++)
			{
				Buffer += " & ";
				if ((l <= Maxl ? CS[k][l] == 0.0 : true)) Buffer += "     ";
				else 
				{
					NB = QString::number(CS[k][l], 'g', (nD[k][l] > 0 ? nD[k][l] + 1 : 8));
					if ((i = NB.indexOf("e")) != -1) 
						Buffer += (lang == 1 ? NB.replace('.', ',') : NB).left(i) + "$\\times 10^{" + 
								QString::number(NB.right(NB.length() - i - 1).toInt()) + "}$";
					else Buffer += (lang == 1 ? NB.replace('.', ',') : NB);
				}
			}
			L << Buffer + "\\\\";
		}
		L << "\\hline";
	}
	/*if (T != 0.0)
	{
		L << "\\hline\\noalign{\\smallskip}";
		L << "$T$ & " + QString::number(T, 'f', (nD[0][0] > 0 ? nD[0][0] - 4 : 9)) + " &     &   \\\\";
		L.append("\\hline");
	}*/
	L << "\\end{tabular}" << "\\end{table}";
	Destroy(DunList, nC);
	Destroy(CS, Maxk + 1);
	Destroy(nD, Maxk + 1);
	if (Err != 0) Destroy(Err, Maxk + 1);
	delete[] RowN;
	return L.join("\n");
}

void DunTable::itemChanged(int row, int /*column*/)
{
	if (BlockIC) return;
	int r, c, NC = Tab->columnCount();
	if (row >= numCoefficients) numCoefficients = row + 1;
	for (r=0; r<= row; r++) for (c=0; c < NC; c++) if (Tab->item(r, c) == 0) 
				Tab->setItem(r, c, new QTableWidgetItem());
	Changed();
}

int DunTable::getMaxk()
{
	int i, B, R=0;
	for (i=0; i<numCoefficients; i++) if ((B = Tab->item(i, 0)->text().toInt()) > R) R = B;
	return R;
}

int DunTable::getMaxl()
{
	int i, B, R=0;
	for (i=0; i<numCoefficients; i++) if ((B = Tab->item(i, 1)->text().toInt()) > R) R = B;
	return R;
}

void DunTable::eIsoO(int *IO, int &NUI, int NI, bool *OI)
{
	int n;
	for (n=0, NUI = 0; n < NI; n++) if (OI[n]) IO[NUI++] = n;
}

double ***DunTable::cvF(int NIso, IsoTab *IsoT, int maxv, int maxk)
{
	double S, T, ***R = Create(NIso, maxv + 1, maxk + 1);
	int I, v, k;
	for (I=0; I < NIso; I++) for (v=0, S = 0.5 * (T = IsoT->rootRRM[I]); v <= maxv; v++, S+=T)
			for (k=1, R[I][v][0] = 1.0; k <= maxk; k++) R[I][v][k] = S * R[I][v][k-1];
	return R;
}

double ***DunTable::cJF(int NIso, IsoTab *IsoT, int maxJ, int maxl, int Omega)
{
	double S, M, ***R = Create(NIso, maxJ + 1, maxl + 1);
	int I, J, l;
	for (I=0; I < NIso; I++) 
		for (J=0, S = - Omega * Omega * (M = IsoT->relRedMass[I]); J <= maxJ; S += 2.0 * (++J) * M)
			for (l=1, R[I][J][0] = 1.0; l <= maxl; l++) R[I][J][l] = S * R[I][J][l-1];
	return R;
}

double ****DunTable::cSpinRotJF(IsoTab *IsoT, int maxN, int maxl)
{
	double S, ****R = Create(2, IsoT->numIso, maxN + 1, maxl + 1);
	int I, N, l;
	for (I=0; I < IsoT->numIso; I++) for (N=0, S = 0.0; N <= maxN; N++, S += IsoT->rootRRM[I])
		for (l=1, R[0][I][N][0] = 1.0, R[1][I][N][0] = -1.0; l <= maxl; l++) 
	{
		R[0][I][N][l] = (l==1? S : S * (S + IsoT->rootRRM[I])) * R[0][I][N][l-1];
		R[1][I][N][l] = (l==1 ? S + IsoT->rootRRM[I] : S * (S + IsoT->rootRRM[I])) * R[1][I][N][l-1];
	}
	return R;
}

void DunTable::getStates(int &N, int &nP, int &eC, Transition **ES, int &gC, Transition **GS, int &mv,
						  int &mJ)
{
	N=0;
	if (molecule == 0 || State == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
								tr("The Dunham coefficient set has to be assigned to a molecule and an electronic state to be fitted!"));
		return;
	}
	int n, m, nT = molecule->getNumTransitions(), *Mv = getmaxv(), MJ = getJMax();
	eC = gC = mv = mJ = 0;
	if (MJ == 0) MJ = cMaxJ;
	for (nP = n = 0; n < nT; n++) 
	{
		if ((ES[eC] = GS[gC] = molecule->getTransitionP(n))->getLowerState() == State)
		{
			if (GS[gC]->getLineTable() != 0) 
			{
				N += GS[gC]->getLineTable()->getNgL(Mv, MJ);
				if ((m = GS[gC]->getLineTable()->getMaxvss()) > mv) mv = m;
				if ((m = GS[gC]->getLineTable()->getMaxJss()) > mJ) mJ = m;
				//GS[gC]->Lines->getObsIso(oIso, NI);
				//for (m=0; m < NI; m++) if (oIso[m]) OIso[m] = true;
				nP += GS[gC++]->getLineTable()->getNgTE();
			}
		}
		else if (ES[eC]->getUpperState() == State) if (ES[eC]->getLineTable() != 0) 
		{
			if ((m = ES[eC]->getLineTable()->getMaxvs()) > mv) mv = m;
			if ((m = ES[eC]->getLineTable()->getMaxJs()) > mJ) mJ = m;
			//ES[eC]->Lines->getObsIso(oIso, NI);
			//for (m=0; m < NI; m++) if (oIso[m]) OIso[m] = true;
			N += ES[eC++]->getLineTable()->getNgTE(Mv, MJ);
		}
	}
	if (N==0) QMessageBox::information(this, tr("QT4MolSpektAn"), 
	  tr("There are no lines or term energies available the Dunham coefficient set could be fitted to!"));
	if (mv != 0) delete[] Mv;
	printf("Ende getStates\n");
}

int **DunTable::getDunList(int &nC, int &NC, bool &FS, int &mk, int &ml)
{
	int n, m;
	QString B;
	if (Tab->rowCount() < numCoefficients) numCoefficients = Tab->rowCount();
	if (ml==0 || mk==0) 
	{
		NC = numCoefficients;
		FS = true;
	}
	else 
	{
		NC = (ml + 1) * (mk + 1);
		FS = false;
	}
	int **DunList = CreateInt(NC, 3);
	if (FS) 
	{
		for (n = nC = 0; n < NC; n++) 
		{
			if ((DunList[nC][0] = Tab->item(n, 0)->text().toInt()) > mk) mk = DunList[nC][0];
			if ((DunList[nC][1] = Tab->item(n, 1)->text().toInt()) > ml) ml = DunList[nC][1];
			B = Tab->item(n, 2)->text();
			DunList[nC][2] = (B == "coeff" || (m = B.toInt()) == 1 ? 0 : (B == "ldcorr" || m == 6 ? 1 : (B == "SpinRot" ? 2 
							: (B == "adCorr" ? 3 : -1))));
			for (m=0; m < nC && DunList[nC][2] != -1; m++) 
				if (DunList[m][0] == DunList[nC][0] && DunList[m][1] == DunList[nC][1] && DunList[m][2] == DunList[nC][2])
					DunList[nC][2] = -1;
			if (DunList[nC][2] != -1) nC++;
		}
		if (nC < NC) NC = numCoefficients = nC;
	}
	else for (n = nC = 0; n < numCoefficients; n++)
	{
		DunList[nC][0] = Tab->item(n, 0)->text().toInt();
		DunList[nC][1] = Tab->item(n, 1)->text().toInt();
		B = Tab->item(n, 2)->text();
		DunList[nC][2] = (B == "coeff" || (m = B.toInt()) == 1 ? 0 : (B == "ldcorr" || m == 6 ? 1 : (B == "SpinRot" ? 2 
							: (B == "adCorr" ? 3 : -1)))); 
		if (DunList[nC][0] <= mk && DunList[nC][1] <= ml && DunList[nC][2] != -1) nC++;
	}
	return DunList;
}

ElState *DunTable::getElState()
{
	return State;
}

void DunTable::Fit(int mk, int ml, bool RemD, int vc1, int Jc1, int vc2, int Jc2)
{
	//printf("Beginn DunTable::Fit\n");
	int n, N, nP, nC, mv, mJ, NI, NC;
	double FQS;
	bool FS;
	if ((vc1 != 0 || vc2 != 0) && Jc1 != Jc2)
	{
		vp1 = vc1;
		Jp1 = Jc1;
		vp2 = vc2;
		Jp2 = Jc2;
	}
	//Transition *GS[MaxStates], *ES[MaxStates];
	FitData *FD = (State != 0 ? State->getFitData() : 0);
	if (FD == 0) 
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"To be able to fit the dunham coefficients, the table has to be assigned to an electronic state which has a fit dataset available!");
		return;
	}
	IsoTab *Iso = molecule->getIso();
	NI = Iso->numIso;
	mJ = FD->getMaxJ();
	mv = FD->getMaxv();
	if ((n = getJMax()) < mJ && n != -1) mJ = n;
	else setJMax(mJ);
	int *Mv = getmaxv();
	if (Mv[0] < mv) mv = Mv[0];
	else if (Mv[0] > mv) setvMax(mv);
	for (n=0; n <= mJ; n++) Mv[n] = mv;
	N = FD->getNumLines(Mv, mJ);
	nP = FD->getNumProgressions(Mv, mJ);
	delete[] Mv;
	//int IsoO[NI];
	//bool oIso[NI], OIso[NI];
	//for (n=0; n < NI; n++) OIso[n] = false; 
	//getStates(N, nP, eC, ES, gC, GS, mv, mJ);
	//if (N==0) return;
	int **DunList = getDunList(nC, NC, FS, mk, ml);
	//eIsoO(IsoO, NI, NI, OIso);
	//printf("ml=%d, mk=%d\n", ml, mk);
	double ***JF = cJF(NI, Iso, mJ, ml, State->getLambda());
	double ***vF = cvF(NI, Iso, mv, mk);
	//printf("NC=%d\n", NC);
	double Res[NC], err[NC];
	int SA[NC], iB;
	//printf("Vor fitDunCoeff\n");
	if (FS) FQS = fitDunCoeff(N, nP, FD, DunList, ml, mk, nC, vF, JF, Res, err, Iso);
	else FQS = fitDunSet(N, nP, FD, DunList, ml, mk, nC, vF, JF, Res, err, RemD, Iso);
	delete Iso;
	//printf("Nach fitDunCoeff\n");	
	for (n=0; n < nC; n++) SA[n] = n;
	for (iB = 0; iB != -1; ) for (n=1, iB = -1; n < nC; n++)
		if (DunList[SA[n-1]][2] > DunList[SA[n]][2] || (DunList[SA[n-1]][2] == DunList[SA[n]][2] 
			&& (DunList[SA[n-1]][1] > DunList[SA[n]][1] || (DunList[SA[n-1]][1] == DunList[SA[n]][1] 
			&& DunList[SA[n-1]][0] > DunList[SA[n]][0])))) 
	{
		iB = SA[n];
		SA[n] = SA[n-1];
		SA[n-1] = iB;
	}
	BlockIC = true;
	if (Tab->columnCount() < 5) Tab->setColumnCount(5);
	for (n=0; n < nC && n < numCoefficients; n++)
	{
		Tab->item(n, 0)->setText(QString::number(DunList[SA[n]][0]));
		Tab->item(n, 1)->setText(QString::number(DunList[SA[n]][1]));
		Tab->item(n, 2)->setText(DunList[SA[n]][2] == 0 ? "coeff" : (DunList[SA[n]][2] == 1 ? "ldcorr" : 
									(DunList[SA[n]][2] == 2 ? "SpinRot" : (DunList[SA[n]][2] == 3 ? "adCorr" : ""))));
		Tab->item(n, 3)->setText(QString::number(Res[SA[n]], 'g', 12));
		if (Tab->item(n, 4) != 0) Tab->item(n, 4)->setText(QString::number(err[SA[n]], 'g', 12));
		else Tab->setItem(n, 4, new QTableWidgetItem(QString::number(err[SA[n]], 'g', 12)));
	}
	Tab->setRowCount(nC);
	for (n = numCoefficients; n < nC; n++)
	{
		Tab->setItem(n, 0, new QTableWidgetItem(QString::number(DunList[SA[n]][0])));
		Tab->setItem(n, 1, new QTableWidgetItem(QString::number(DunList[SA[n]][1])));
		Tab->setItem(n, 2, new QTableWidgetItem(DunList[SA[n]][2] == 0 ? "coeff" : (DunList[SA[n]][2] == 1 ? "ldcorr" : 
									(DunList[SA[n]][2] == 2 ? "SpinRot" : ""))));
		Tab->setItem(n, 3, new QTableWidgetItem(QString::number(Res[SA[n]], 'g', 12)));
		Tab->setItem(n, 4, new QTableWidgetItem(QString::number(err[SA[n]], 'g', 12)));
	}
	numCoefficients = nC;
	QMessageBox::information(this, "MolSpektAnalysis", 
							"The resulting sigma is " + QString::number(sqrt(FQS / double(N - nC)), 'f', 2));
	Destroy(DunList, NC);
	BlockIC = false;
	Changed();
}

void DunTable::calcXY(Transition **GS, int gC, Transition **ES, int eC, int NP, int xx[][4], 
					  double *yy, double *ssig, double *SS)
{	
	TableLine *TL;
	TermEnergy *TE;
	int n, l, c, m, *mv = getmaxv(), mJ = getJMax();
	if (mJ == 0) mJ = cMaxJ;
	//printf("eC=%d, gC=%d\n", eC, gC);
	for (n=0; n < NP; n++) SS[n] = 0.0;
	for (n = l = 0; n < gC; n++)
	{
		GS[n]->getLineTable()->SortProg();
		GS[n]->getLineTable()->getgoodLines(c, TL, mv, mJ);
		if (c == 0) continue;
		/*if (ES[n]->upperState == 0) continue;
		if (ES[n]->upperState->termTable == 0) continue;
		UData = ES[n]->upperState->termTable->getData();
		uNIso = ES[n]->upperState->termTable->getNumIso();
		uMv = ES[n]->upperState->termTable->getMaxv();
		uMJ = ES[n]->upperState->termTable->getMaxJ();*/
		for (m=0; m<c; m++)
		{
			xx[l][0] = TL[m].Iso;
			xx[l][1] = TL[m].vss;
			xx[l][2] = TL[m].Jss;
			xx[l][3] = TL[m].PN - 1;
			yy[l] = TL[m].WN;
			SS[TL[m].PN - 1] += 1.0 / (ssig[l++] = TL[m].err * TL[m].err);
		}
		delete[] TL;
	}
	//printf("Nach gC\n");
	for (n=0; n < eC; n++)
	{
		ES[n]->getLineTable()->getgoodTE(c, TE, mv, mJ);
		//printf("Nach getgoodTE, c=%d\n", c);
		if (c == 0) continue;
		for (m=0; m<c; m++)
		{
			xx[l][0] = TE[m].Iso;
			xx[l][1] = TE[m].v;
			xx[l][2] = TE[m].J;
			xx[l][3] = -1;
			yy[l] = TE[m].E;
			//printf("xx[%d][0]=%d, xx[%d][1]=%d, xx[%d][2]=%d, xx[%d][3]= %d, y[%d]=%f\n", 
				//   l, xx[l][0], l, xx[l][1], l, xx[l][2], l, xx[l][3], l, yy[l]);
			ssig[l++] = 0.0001;
		}
		delete[] TE;
	}
	//printf("Ende von calcXY\n");
	if (mv != 0) delete[] mv;
}

void DunTable::calcXY(FitData* FD, int NP, int &N, int xx[][5], double* yy, double* ssig, 
					  int *RN, double* SS)
{
	TableLine *TL;
	int n, l, c, m, *mv = getmaxv(), mJ = getJMax(), *RNi;
	if (mJ == 0) mJ = cMaxJ;
	for (n=0; n < NP; n++) SS[n] = 0.0;
	FD->getData(TL, c, RNi, sortByProg, mv, mJ);
	for (m=l=0, n=-1; m<c && l<N; m++)
	{
		xx[l][0] = TL[m].Iso;
		xx[l][1] = TL[m].vss;
		xx[l][2] = TL[m].Jss;
		xx[l][4] = 2 * TL[m].FC + (TL[m].Js == TL[m].Jss ? 1 : 0);
		yy[l] = TL[m].WN;
		ssig[l] = TL[m].err * TL[m].err;
		RN[l] = RNi[m];
		if (TL[m].isTE) xx[l++][3] = -1;
		else
		{
			if (n > -1 ? TL[m].PN != TL[m-1].PN || TL[m].LTab != TL[m-1].LTab : true)
			{
				if (TL[m].PN == TL[m+1].PN) n++;
				else continue;
				printf("%d, ", m);
			}
			xx[l][3] = n;
			SS[n] += 1.0 / ssig[l++];
		}
	}
	printf("\n");
	N=l;
	if (c>0)
	{
		delete[] TL;
		delete[] RNi;
	}
	if (mv != 0) delete[] mv;
}

double DunTable::fitDunCoeff(int N, int nP, FitData *FD, int **DunList, int /*Ml*/, 
							 int /*Mk*/, int nC, double ***vF, double ***JF, 
							 double *Res, double *err, IsoTab *Iso)
{
	//printf("nP=%d, N=%d\n", nP, N);	
    double *SS = new double[nP], *ssig = new double[N], *yy = new double[N], *dev = new double[N];
    int (*xx)[5] = new int[N][5], *RowN = new int[N];
	double Ret;
	calcXY(FD, nP, N, xx, yy, ssig, RowN, SS);
	//printf("Nach calcXY\n");
	if (N==0) return 0.0;
	Ret = fitDun(N, nC, DunList, vF, JF, Res, err, xx, yy, ssig, SS, dev, Iso);
	FD->setDev(dev, RowN, N);
    delete[] SS;
    delete[] ssig;
    delete[] yy;
    delete[] dev;
    delete[] xx;
    delete[] RowN;
	return Ret;
}

double DunTable::fitDunSet(int N, int nP, FitData *FD, int **DunList, int Ml, int Mk, 
						   int &nC, double ***vF, double ***JF, 
		                   double *Res, double *err, bool RemD, IsoTab *Iso)
{
	double SS[nP], ssig[N], yy[N], R, r, aFQS = 0.0, FQS = 0.0, dev[N], md;
	bool success = false;
	int xx[N][5], n, m, k=Mk+1, l=Ml, M=0, p=nC, s=-1, count=0, b=-1, RN[N];
	calcXY(FD, nP, N, xx, yy, ssig, RN, SS);
	int **D = CreateInt(Mk+1, Ml+1);
	for (n=0; n <= Ml; n++) for (m=0; m <= Mk; m++) D[m][n] = 0;
	for (n=0; n < nC; n++) D[DunList[n][0]][DunList[n][1]] = 1;
	printf("Vor main loop, nC=%d\n", nC);
	//RemD = true;
	while (true)
	{
		printf("M=%d, nC=%d, k=%d, l=%d, p=%d, s=%d\n", M, nC, k, l, p, s);
		if (M==1)
		{
			if ((FQS < 1.01 * aFQS || R < 10.0) && nC > 0 && !isnan(FQS) && !isnan(R))
			{
				if (FQS != aFQS) printf("Successfull!\n");
				for (n=0, R=1e99; n<nC; n++) if ((r = fabs(Res[n] / err[n])) < R) 
				{
					printf("n=%d, r=%f\n", n, r);
					R=r;
					k = DunList[n][0];
					l = DunList[n][1];
				}
				if (R < 1e3 && FQS < N) D[k][l] = -1;
				else M++;
				//printf("vor p loop\n");
				for (p=n=0; n<=l; n++) for (m=0; m <= (n<l ? Mk : k); m++) if (D[m][n] == 1) p++;
				printf("Trying to remove coefficient (k=%d,l=%d) with r=%f\n", k, l, R);
			}
			else 
			{
				if (FQS > aFQS) 
				{
					printf("Not successfull!\n");
					D[k][l]=1;
					FQS = aFQS;
				}
				M++;
			}
		}
		if (M==3)
		{
			r = fabs(Res[p] / err[p]);
			if (r < 10.0 || (FQS > 0.99 * aFQS && r < 1e3) || isnan(r) || isnan(FQS)) 
			{
				printf("Not successful, r=%f!\n", r);
				D[k][l] = -1;
				s=p;
				FQS = aFQS;
			}
			else 
			{
				printf("Successful, r=%f!\n", r);
				M=1;
				aFQS = FQS;
				success = true;
				//for (k=0; k<=Mk; k++) for (l=0; l<=Ml; l++) if (D[k][l]==-1) D[k][l]=0;
				continue;
			}
		}
		if (M==2)
		{
			printf("Mode=2\n");
			if (nC > 0) b=p;
			k=l=p=0;
			//if (eC == 0) k=1;
			M++;
		}
		if (M==3)
		{
			count = 0;
			while (D[k][l] != 0) 
			{
				printf("D[%d][%d]=%d\n", k, l, D[k][l]);
				if (D[k][l]==1) p++;
				if (++k > Mk)
				{
					if (l < Ml)
					{
						l++;
						k=0;
					}
					else 
					{
						if (s!=-1) 
						{
							M=0;
							break;
						}
						if (RemD)
						{
							for (k=0, md = 0.0; k<N; k++) 
							{
								if (dev[k] > md) md = dev[l=k];
								if (ssig[k] > 100.0) if(ssig[k] * dev[k] / (ssig[k] - 100.0) < 9.0) 
								{
									SS[xx[k][3]] += 100.0 / (ssig[k] * (-100.0 + ssig[k]));
									ssig[k] -= 100.0;
									printf("A:setting ssig[%d]=%f for Iso=%d, v=%d, J=%d\n", 
										   k, ssig[k], xx[k][0], xx[k][1], xx[k][2]);
								}
							}
							if (md > 9.0) 
							{
								SS[xx[l][3]] -= 100.0 / (ssig[l] * (100.0 + ssig[l]));
								ssig[l] += 100.0;
								printf("R:setting ssig[%d]=%f for Iso=%d, v=%d, J=%d\n", 
									   	l, ssig[l], xx[l][0], xx[l][1], xx[l][2]);
								success = true;
								M=0;
								aFQS = FQS;
								break;
							}
						}
						if (success)
						{
							M=1;
							for (k=0; k<=Mk; k++) for (l=0; l<=Ml; l++) if (D[k][l]==-1) D[k][l]=0;
							success = false;
							aFQS = FQS;
							break;
						}
						else 
						{
							M = 4;
							k=l=p=0;
							break;
						}
					}
				}
			}
			if (k > Mk) continue;
			if (M<4)
			{
				D[k][l] = 1;
				printf("Adding coefficient (k=%d, l=%d)\n", k, l);
			}
		}
		if (b>0)
		{
			s=b-1;
			b=-1;
		}
		if (s==-1) 
		{
			s=p;
			m=k;
			n=l;
		}
		else
		{
			m = DunList[s][0];
			n = DunList[s][1];
		}
		if (b==0)
		{
			b=-1;
			n=m=s=0;
		}
		for (nC = s; m<=Mk || n<Ml; m++)
		{
			if (m>Mk)
			{
				n++;
				m=0;
			}
			if (D[m][n] == 1)
			{
				DunList[nC][0] = m;
				DunList[nC++][1] = n;
			}
		}
		//if ((DunList[0][0] != 0 || DunList[0][1] != 0) && nC > 0) break;
		if (aFQS != FQS)
		{
			count=0;
			aFQS = FQS;
		}
		for (n=0; n<nC; n++) printf("DunList[%d]=(%d,%d)\n", n, DunList[n][0], DunList[n][1]);
		FQS = fitDun(N, nC, DunList, vF, JF, Res, err, xx, yy, ssig, SS, dev, Iso);
		printf("FQS=%g\n", FQS);
		if (M==4) break;
		if (aFQS == 0.0) aFQS = FQS;
		if (M==0) 
		{
			M++;
			aFQS = FQS;
		}
		s=-1;
		if (count++ == 5) 
		{
			printf("Error: endless loop detected!\n");
			break;
		}
	}
	Destroy(D, Mk+1);
	QFile FR("Rem.log"), FU("Res.log");
	FR.open(QIODevice::WriteOnly);
	FU.open(QIODevice::WriteOnly);
	QTextStream SR(&FR), SU(&FU);
	SR << "PN	Iso	v	J	obs-calc\n";
	SU << "PN	Iso	v	J	obs-calc\n";
	for (k=0; k<N; k++)
	{
		if (ssig[k] < 100.0) SU << xx[k][3] << "\t" << xx[k][0] << "\t" << xx[k][1] 
	      						<< "\t" << xx[k][2] << "\t" 
								<< QString::number(sqrt(ssig[k] * dev[k]), 'f', 4) << "\n";
		else SR << xx[k][3] << "\t" << xx[k][0] << "\t" << xx[k][1] 
				<< "\t" << xx[k][2] << "\t" 
				<< QString::number(sqrt(ssig[k] * dev[k]), 'f', 4) << "\n";
	}
	FR.close();
	FU.close();
	printf("Finished!\n");
	return FQS;
}

void DunTable::calcRe(double& Re, double& Err)
{
    double Be;
    int n = getBe(Be);
	IsoTab *Iso = (molecule != 0 ? molecule->getIso() : 0);
    if (Iso == 0 || Be <= 0.0)
	{
		Re = Err = 0.0;
		return;
	}
    if (Be > 0.0)
    {
        Re = 1e9 / M_PI * sqrt(C_h / (8.0 * Iso->redMass[Iso->refIso] * Be * C_u * C_c));
        Err = 0.5 * Re / Be * Tab->item(n, 4)->text().toDouble();
    }
}

int DunTable::getBe(double &o_Be) const
{
    int n;
    QString Type;
    for (n=0; n < numCoefficients; n++)
    {
        Type = Tab->item(n, 2)->text();
        if (Tab->item(n, 0)->text().toInt() == 0 && Tab->item(n, 1)->text().toInt() == 1 && (Type == "coeff" || Type.toInt() == 1))
        {
            o_Be = Tab->item(n, 3)->text().toDouble();
            return n;
        }
    }
    o_Be = 0.0;
    return 0;
}

void DunTable::testKratzer(double& Y02f, double& Y02c)
{
	double Be = 0.0, we = 0.0;
	int n, k, l, t;
	Y02f = 0.0;
	Y02c = 0.0;
	QString Type;
	for (n=0; n < numCoefficients; n++) 
	{
		Type = Tab->item(n, 2)->text();
		if (Type == "coeff" || (t = Type.toInt()) == 1)
		{
			k = Tab->item(n, 0)->text().toInt();
			l = Tab->item(n, 1)->text().toInt();
			if (k == 1 && l == 0) we = Tab->item(n, 3)->text().toDouble();
			else if (k == 0 && l == 1) Be = Tab->item(n, 3)->text().toDouble();
			else if (k == 0 && l == 2) Y02f = Tab->item(n, 3)->text().toDouble();
		}
	}
	if (we > 0.0 && Be > 0.0) Y02c = -4.0 * Be * Be * Be / (we * we);
}

void DunTable::calcTeY00(double& Te, double& TeErr, double& Y00, double& Y00Err)
{
	double TeY00 = 0.0, Y20 = 0.0, Y11 = 0.0, Be = 0.0, we = 0.0, TeY00Err = -1.0, Y20Err = -1.0, Y11Err = -1.0, BeErr = -1.0, weErr = -1.0;
	double d6Beq, a1, a2, a1Err, a2Err, T1, T2, T3, d3Be;
	int n, k, l;
	QString Type;
	for (n=0; n < numCoefficients; n++)
	{
		Type = Tab->item(n, 2)->text();
		if (Type == "coeff" || Type.toInt() == 1)
		{
			if ((k = Tab->item(n, 0)->text().toInt()) == 0)
			{
				if ((l = Tab->item(n, 1)->text().toInt()) == 1)
				{
					Be = Tab->item(n, 3)->text().toDouble();
					BeErr = Tab->item(n, 4)->text().toDouble();
				}
				else if (l==0)
				{
					TeY00 = Tab->item(n, 3)->text().toDouble();
					TeY00Err = Tab->item(n, 4)->text().toDouble();
				}
			}
			else if (k==1)
			{
				if ((l = Tab->item(n, 1)->text().toInt()) == 0)
				{
					we = Tab->item(n, 3)->text().toDouble();
					weErr = Tab->item(n, 4)->text().toDouble();
				}
				else if (l==1)
				{
					Y11 = Tab->item(n, 3)->text().toDouble();
					Y11Err = Tab->item(n, 4)->text().toDouble();
				}
			}
			else if (k==2 && Tab->item(n, 1)->text().toInt() == 0)
			{
				Y20 = Tab->item(n, 3)->text().toDouble();
				Y20Err = Tab->item(n, 4)->text().toDouble();
			}
		}
	}
	if (Y20 == 0.0 || Be == 0.0 || we == 0.0 || Y11 == 0.0)
	{
		Te = TeErr = Y00 = Y00Err = 0.0;
		return;
	}
	d3Be = 1.0 / (3.0 * Be);
	d6Beq = 1.5 * d3Be * d3Be;
	a1 = Y11 * we * d6Beq - 1.0;
	T1 = Y11Err * we;
	T2 = Y11 * weErr;
	T3 = Y11 * we * BeErr * 6.0 * d3Be * d6Beq;
	a1Err = d6Beq * sqrt(T1 * T1 + T2 * T2 + T3 * T3);
	a2 = 2.0 * Y20 * d3Be + 1.25 * a1 * a1;
	T1 = 2.0 * d3Be * Y20Err;
	T2 = 4.0 * Y20 * BeErr * d6Beq;
	T3 = 2.5 * a1 * a1Err;
	a2Err = sqrt(T1 * T1 + T2 * T2 + T3 * T3);
	Y00 = Be * (0.375 * a2 - 0.21875 * a1 * a1);
	T1 = 0.4375 * a1 * a1Err;
	Y00Err = Be * sqrt(0.140625 * a2Err * a2Err + T1 * T1);
	Te = TeY00 - Y00;
	TeErr = sqrt(TeY00Err * TeY00Err + Y00Err * Y00Err);
}
