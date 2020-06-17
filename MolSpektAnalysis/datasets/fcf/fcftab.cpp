//
// C++ Implementation: FCFTab
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "fcftab.h"
#include "MainWindow.h"
#include "molecule.h"
#include "potential.h"
#include "utils.h"
#include "isotab.h"
#include "CoupledSineWaveFunc.h"
#include "Transition.h"

#include <QComboBox>
#include <QDateTime>
#include <QFile>

#include <math.h>


FCFTab::FCFTab(MainWindow* MW, Molecule* M): TableWindow(FranckCondonTable, MW, M), m_NChannelsQ(0), m_NChannelsPR(0), m_ChannelFQ(NULL), m_ChannelFPR(NULL)
{
	NIso = 0;
	NJ = 0;
	Nvs = 0;
	Nvss = 0;
	Data = 0;
	aPQR = 1;
	Block = false;
	St1 = St2 = 0;
	Trans = 0;
	IsoTab *IsoT = M->getIso();
	int i;
	for (i=0; i < (nIso = IsoT->numIso); i++) Iso->addItem(IsoT->getIsoName(i));
	connect(Iso, SIGNAL(currentIndexChanged(int)), this, SLOT(IsoChanged(int)));
	connect(JsB, SIGNAL(currentIndexChanged(int)), this, SLOT(JsChanged()));
	connect(JssB, SIGNAL(currentIndexChanged(int)), this, SLOT(JssChanged()));
	Iso->setCurrentIndex(IsoT->refIso);
	delete IsoT;
}

FCFTab::~FCFTab()
{
	DestroyData();
}

bool FCFTab::calcFCF(ElState* State1, ElState* State2, int Maxvs, int Maxvss, int MaxJ, int NumWFPoints)
{
    Potential *pot1 = State1->getPotential(), *pot2 = State2->getPotential();
    if (pot1 == 0 || pot2 == 0) return false;
    CoupledSineWaveFunc* SWF = pot1->getCoupledSineWaveFunc();
    Molecule* mol = State1->getMolecule();
    double *ChannelFQ, *ChannelFPR;
    int NChannelsQ = 0, NChannelsPR = 0;
    if (0 != SWF)
    {
        int NC = SWF->getNChannels();
        ChannelFQ = new double[NC];
        ChannelFPR = new double[NC];
        ElState* lg1State = nullptr;
        for (int n = 0; n < NC; ++n)
        {
            int c;
            ElState* St;
            SWF->getStateComp(n, St, c);
            if (c == 0)
            {
                ChannelFPR[NChannelsPR++] = mol->getTransitionStrength(St, State2);
                if (St->getLambda() >= 1) lg1State = St;
            }
            else ChannelFQ[NChannelsQ++] = mol->getTransitionStrength(St, State2);
        }
        if (0 == NChannelsQ && lg1State != nullptr)
        {
            NChannelsQ = 1;
            ChannelFQ[0] = mol->getTransitionStrength(lg1State, State2);
        }
    }
    else
    {
        if (State1->getLambda() >= 1)
        {
            ChannelFQ = new double[++NChannelsQ];
            ChannelFQ[0] = mol->getTransitionStrength(State1, State2);
        }
        ChannelFPR = new double[++NChannelsPR];
        ChannelFPR[0] = mol->getTransitionStrength(State1, State2);
    }
    m_ChannelFQ = ChannelFQ;
    m_ChannelFPR = ChannelFPR;
    m_NChannelsQ = NChannelsQ;
    m_NChannelsPR = NChannelsPR;
	int Nv1, Nv2, Nv3, Nv4, *P1 = 0, *P2 = 0, *P3 = 0, *P4 = 0, *ri1 = 0, *ri2 = 0, *ri3 = 0;
	int *ri4 = 0, *ra1 = 0, *ra2 = 0, *ra3 = 0, *ra4 = 0;
	int Iso, i, JStart1 = State1->getJStart(0, 0), JStart2 = State2->getJStart(0, 0);
	int JStep, nJ = MaxJ + 1, J, Js;
    double *E1 = 0, *E2 = 0, *E3 = 0, *E4 = 0, **WF2 = 0;
	double *F1 = 0, *F2 = 0, *F3 = 0, *F4 = 0;
    double *N1 = 0, *N2 = 0, *N3 = 0, *N4 = 0;
	Pot1->setText(pot1->getName());
	Pot2->setText(pot2->getName());
	Date->setText(QDateTime::currentDateTime().toString("dd.MM.yyyy hh:mm"));
	St1 = State1;
	St2 = State2;
	NIso = new int[3];
	NJ = new int*[3];
	Nvss = new int**[3];
	Nvs = new int**[3];
	Data = new float****[3];
    if (State1->getLambda() == 0 && SWF == nullptr)
	{
        double **WF1 = 0, ***FCF = 0;
		NIso[1] = 0;
		NIso[0] = NIso[2] = nIso;
		NJ[0] = new int[nIso];
		NJ[1] = 0;
		NJ[2] = new int[nIso];
		Nvs[0] = new int*[nIso];
		Nvs[2] = new int*[nIso];
		Nvss[0] = new int*[nIso];
		Nvss[2] = new int*[nIso];
		Nvs[1] = Nvss[1] = 0;
		Data[0] = new float***[nIso];
		Data[1] = 0;
		Data[2] = new float***[nIso];
		for (Iso = 0; Iso < nIso; Iso++)
		{
			Nvs[0][Iso] = new int[nJ];
			Nvs[2][Iso] = new int[nJ];
			Nvss[0][Iso] = new int[nJ];
			Nvss[2][Iso] = new int[nJ];
			Data[0][Iso] = new float**[nJ];
			Data[2][Iso] = new float**[nJ];
			for (J=0; J < nJ; J++)
			{
				Nvs[0][Iso][J] = Nvs[2][Iso][J] = Nvss[0][Iso][J] = Nvss[2][Iso][J] = 0;
				Data[0][Iso][J] = Data[2][Iso][J] = 0;
			}
			for (J = Js = JStart2, JStep = molecule->getJStep(Iso);
					true; J+=2)
			{
				if (J > MaxJ || (J > Js && Nv1 == 0))
				{
					if (Js != JStart2 || JStep == 2) break;
					else J = ++Js;
					if (JStart1 < Js)
					{
						cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, Nv1);
                        pot1->getWaveFunc(Iso, J-1, Nv1 = Maxvs, E1, WF1, ri1, ra1, P1, F1, N1, NumWFPoints);
					}
				}
				cleanUp(E2, WF2, ri2, ra2, P2, F2, N2, Nv2);
                pot2->getWaveFunc(Iso, J, Nv2 = Maxvss, E2, WF2, ri2, ra2, P2, F2, N2, NumWFPoints);
				if (J > JStart1 && Nv2 > 0)
				{
					pot1->getFCF(Nv1, Nv2, E1, E2, WF1, WF2, ri1, ri2, ra1, ra2, P1, P2,
								 F1, F2, N1, N2, FCF);
					setData(0, Iso, J, Nv2, Nv1, FCF);
				}
				if (J+1 > MaxJ || Nv2 == 0)
				{
					if (Js != JStart2 || JStep == 2) break;
					else J = ++Js;
				}
				cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, Nv1);
                pot1->getWaveFunc(Iso, J+1, Nv1 = Maxvs, E1, WF1, ri1, ra1, P1, F1, N1, NumWFPoints);
				if (Nv1 > 0)
				{
					pot1->getFCF(Nv1, Nv2, E1, E2, WF1, WF2, ri1, ri2, ra1, ra2, P1, P2,
								 F1, F2, N1, N2, FCF);
					setData(2, Iso, J, Nv2, Nv1, FCF);
				}
			}
			for (NJ[0][Iso] = J; Data[0][Iso][NJ[0][Iso] - 1] == 0; NJ[0][Iso]--) ;
			for (NJ[2][Iso] = J; Data[2][Iso][NJ[2][Iso] - 1] == 0; NJ[2][Iso]--) ;
		}
        cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, Nv1);
	}
	else
	{
        double ***WF1 = 0, ***WF3 = 0, ***WF4 = 0, ****FCF = 0;
        int c;
		for (i=0; i<3; i++)
		{
			NIso[i] = nIso;
			NJ[i] = new int[nIso];
			Nvs[i] = new int*[nIso];
			Nvss[i] = new int*[nIso];
			Data[i] = new float***[nIso];
		}
		for (Iso = 0; Iso < nIso; Iso++)
		{
			for (i=0; i<3; i++)
			{
				Nvs[i][Iso] = new int[nJ];
				Nvss[i][Iso] = new int[nJ];
				Data[i][Iso] = new float**[nJ];
				for (J=0; J < nJ; J++)
				{
					Nvs[i][Iso][J] = Nvss[i][Iso][J] = 0;
					Data[i][Iso][J] = 0;
				}
			}
			if ((JStep = molecule->getJStep(Iso)) == 2) for (J = JStart2; J <= MaxJ; J += JStep)
			{
                cleanUp(E2, WF2, ri2, ra2, P2, F2, N2, Nv2);
                pot2->getWaveFunc(Iso, J, Nv2 = Maxvss, E2, WF2, ri2, ra2, P2, F2, N2, NumWFPoints);
				if (Nv2 == 0) break;
				if (J > JStart1)
				{
                    for (c=0, FCF = new double***[m_NChannelsPR]; c < m_NChannelsPR; ++c)
                        pot1->getFCF(Nv1, Nv2, E1, E2, WF1[c], WF2, ri1, ri2, ra1, ra2, P1, P2, F1, F2, N1, N2, FCF[c]);
					setData(0, Iso, J, Nv2, Nv1, FCF);
				}
				if (J >= JStart1)
				{
                    cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, m_NChannelsPR, Nv1);
                    if (SWF == 0 || m_NChannelsQ == 1)
                        pot1->getWaveFunc(Iso, J, Nv1 = Maxvs, E1, (WF1 = new double**[1])[0], ri1, ra1, P1, F1, N1, NumWFPoints);
                    else SWF->getWaveFuncs(Iso, J, Nv1 = Maxvs, m_NChannelsQ, rmin, rmax, NumWFPoints, WF1, ri1, ra1, P1, F1, N1);
					if (Nv1 == 0) break;
                    for (c=0, FCF = new double***[m_NChannelsQ]; c < m_NChannelsQ; ++c)
                        pot1->getFCF(Nv1, Nv2, E1, E2, WF1[c], WF2, ri1, ri2, ra1, ra2, P1, P2, F1, F2, N1, N2, FCF[c]);
					setData(1, Iso, J, Nv2, Nv1, FCF);
				}
                cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, m_NChannelsQ, Nv1);
                if (SWF == 0 || m_NChannelsPR == 1)
                    pot1->getWaveFunc(Iso, J+1, Nv1 = Maxvs, E1, (WF1 = new double**[1])[0], ri1, ra1, P1, F1, N1, NumWFPoints);
                else SWF->getWaveFuncs(Iso, J+1, Nv1 = Maxvs, m_NChannelsPR, rmin, rmax, NumWFPoints, WF1, ri1, ra1, P1, F1, N1);
				if (Nv1 == 0) break;
                for (c=0, FCF = new double***[m_NChannelsPR]; c < m_NChannelsPR; ++c)
                    pot1->getFCF(Nv1, Nv2, E1, E2, WF1[c], WF2, ri1, ri2, ra1, ra2, P1, P2, F1, F2, N1, N2, FCF[c]);
				setData(2, Iso, J, Nv2, Nv1, FCF);
			}
			else 
			{
				for (J = JStart2; J <= MaxJ; J++)
				{
					cleanUp(E2, WF2, ri2, ra2, P2, F2, N2, Nv2);
                    pot2->getWaveFunc(Iso, J, Nv2 = Maxvss, E2, WF2, ri2, ra2, P2, F2, N2, NumWFPoints);
					if (Nv2 == 0) break;
                    cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, m_NChannelsPR, Nv1);
                    shiftPar(E1, E4, WF1, WF4, ri1, ri4, ra1, ra4, P1, P4, F1, F4, N1, N4, Nv1, Nv4);
					if (J > JStart1)
					{
                        for (c=0, FCF = new double***[m_NChannelsPR]; c < m_NChannelsPR; ++c)
                            pot1->getFCF(Nv1, Nv2, E1, E2, WF1[c], WF2, ri1, ri2, ra1, ra2, P1, P2, F1, F2, N1, N2, FCF[c]);
						setData(0, Iso, J, Nv2, Nv1, FCF);
					}
					if (J >= JStart1)
					{
                        cleanUp(E3, WF3, ri3, ra3, P3, F3, N3, m_NChannelsQ, Nv3);
                        if (SWF == 0 || m_NChannelsQ == 1)
                            pot1->getWaveFunc(Iso, J, Nv3 = Maxvs, E3, (WF3 = new double**[1])[0], ri3, ra3, P3, F3, N3, NumWFPoints);
                        else SWF->getWaveFuncs(Iso, J, Nv3 = Maxvs, m_NChannelsQ, rmin, rmax, NumWFPoints, WF3, ri3, ra3, P3, F3, N3);
                        for (c=0, FCF = new double***[m_NChannelsQ]; c < m_NChannelsQ; ++c)
                            pot1->getFCF(Nv3, Nv2, E3, E2, WF3[c], WF2, ri3, ri2, ra3, ra2, P3, P2, F3, F2, N3, N2, FCF[c]);
						setData(1, Iso, J, Nv2, Nv3, FCF);
					}
                    cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, m_NChannelsPR, Nv1);
                    if (SWF == 0 || m_NChannelsPR == 1)
                        pot1->getWaveFunc(Iso, J+1, Nv4 = Maxvs, E4, (WF4 = new double**[1])[0], ri4, ra4, P4, F4, N4, NumWFPoints);
                    else SWF->getWaveFuncs(Iso, J+1, Nv4 = Maxvs, m_NChannelsPR, rmin, rmax, NumWFPoints, WF4, ri4, ra4, P4, F4, N4);
					if (Nv4 == 0) break;
                    for (c=0, FCF = new double***[m_NChannelsPR]; c < m_NChannelsPR; ++c)
                        pot1->getFCF(Nv4, Nv2, E4, E2, WF4[c], WF2, ri4, ri2, ra4, ra2, P4, P2, F4, F2, N4, N2, FCF[c]);
					setData(2, Iso, J, Nv2, Nv4, FCF);
				}
                cleanUp(E3, WF3, ri3, ra3, P3, F3, N3, m_NChannelsQ, Nv3);
                cleanUp(E4, WF4, ri4, ra4, P4, F4, N4, m_NChannelsPR, Nv4);
			}
			for (i=0; i<3; i++) 
				for (NJ[i][Iso] = J; Data[i][Iso][NJ[i][Iso] - 1] == 0; NJ[i][Iso]--) ;
		}
        cleanUp(E1, WF1, ri1, ra1, P1, F1, N1, m_NChannelsPR, Nv1);
        cleanUp(E3, WF3, ri3, ra3, P3, F3, N3, m_NChannelsQ, Nv3);
	}
	cleanUp(E2, WF2, ri2, ra2, P2, F2, N2, Nv2);
	Changed();
	return true;
}

void FCFTab::cleanUp(double* E, double** WF, int* ri, int* ra, int* P, double* F, double* N, int Nv)
{
	if (E == 0) return; 
	delete[] E;
	E=0;
	Destroy(WF, Nv);
	delete[] ri;
	delete[] ra;
	delete[] P;
	delete[] F;
	delete[] N;
}

void FCFTab::cleanUp(double *E, double ***WF, int *ri, int *ra, int *P, double *F, double *N, const int NChannels, const int Nv)
{
    if (E == 0) return;
    delete[] E;
    E=0;
    Destroy(WF, NChannels, Nv);
    delete[] ri;
    delete[] ra;
    delete[] P;
    delete[] F;
    delete[] N;
}

void FCFTab::DestroyData()
{
	int i, I, J, v;
	if (Data != 0)
	{
		for (i=0; i<3; i++) if (Data[i] != 0)
		{
			for (I=0; I < NIso[i]; I++) if (Data[i][I] != 0)
			{
				for (J=0; J < NJ[i][I]; J++) if (Data[i][I][J] != 0)
				{
					for (v = 0; v < Nvs[i][I][J]; v++) delete[] Data[i][I][J][v];
					delete[] Data[i][I][J];
				}
				delete[] Data[i][I];
				delete[] Nvss[i][I];
				delete[] Nvs[i][I];
			}
			delete[] Data[i];
			delete[] Nvss[i];
			delete[] Nvs[i];
			delete[] NJ[i];
		}
		delete[] Data;
		delete[] Nvss;
		delete[] Nvs;
		delete[] NJ;
		delete[] NIso;
	}
    if (m_ChannelFQ != NULL) delete[] m_ChannelFQ;
    if (m_ChannelFPR != NULL) delete[] m_ChannelFPR;
}

void FCFTab::getFCF(int PQR, int Iso, int Js, int vs, int& N, float*& D)
{
	D=0;
	N=0;
	if (Data != 0) if (Iso < NIso[PQR]) if (Js < NJ[PQR][Iso]) if (vs < Nvs[PQR][Iso][Js])
	{
		N = Nvss[PQR][Iso][Js];
		D = Data[PQR][Iso][Js][vs];
	}
}

void FCFTab::getFCF(int JD, int Iso, int Js, int vs, int& N, double* D)
{
	int n, PQR = (JD == 0 ? 1 : 0), M = N;
	D=0;
	N=0;
	if (Data != 0) if (Iso < NIso[PQR]) if (Js < NJ[PQR][Iso]) if (vs < Nvs[PQR][Iso][Js])
	{
		if ((N = Nvss[PQR][Iso][Js]) > M) N=M;
		if (PQR == 1) for (n=0; n<N; n++) D[n] = Data[PQR][Iso][Js][vs][n];
		else
		{
			for (n=0; n<N; n++) D[n] = (n < Nvss[2][Iso][Js] ? Data[2][Iso][Js][vs][n] : 0.0);
			for (n=0; n<N; n++) D[n+N] = Data[PQR][Iso][Js][vs][n];
		}
	}
}

ElState* FCFTab::getState1()
{
	return St1;
}

ElState* FCFTab::getState2()
{
	return St2;
}

Transition* FCFTab::getTransition()
{
	return Trans;
}

void FCFTab::IsoChanged(int NIso)
{
	int J, n=0, aJ = JsB->currentText().toInt(), d = 1000, bn = 0;
	JsB->clear();
	if (Data != 0) for (J=0; J < NJ[2][NIso]; J++) 
		if ((J>0 ? Nvs[2][NIso][J-1] > 0 : false) || (Nvs[1] != 0 ? Nvs[1][NIso][J] > 0 : false) 
			|| (J+1 < NJ[0][NIso] ? Nvs[0][NIso][J+1] > 0 : false)) 
	{
		JsB->addItem(QString::number(J));
		if (fabs(J - aJ) < d)
		{
			d = fabs(J - aJ);
			bn = n;
		}
		n++;
	}
	JsB->setCurrentIndex(bn);
}

void FCFTab::JsChanged()
{
	int J = JsB->currentText().toInt(), NIso = Iso->currentIndex();
	Block = true;
	JssB->clear();
	if (J>0 ? Nvs[2][NIso][J-1] > 0 : false) JssB->addItem(QString::number(J-1));
	if (Nvs[1] != 0 ? Nvs[1][NIso][J] > 0 : false)
	{
		JssB->addItem(QString::number(J));
		if (aPQR == 1) JssB->setCurrentIndex(JssB->count() - 1);
	}
	if (J+1 < NJ[0][NIso] ? Nvs[0][NIso][J+1] > 0 : false) JssB->addItem(QString::number(J+1));
	if (aPQR == 0) JssB->setCurrentIndex(0);
	else if (aPQR == 2) JssB->setCurrentIndex(JssB->count() - 1);
	Block = false;
	JssChanged();
}

void FCFTab::JssChanged()
{
	if (Block) return;
	int vs, vss, Js = JsB->currentText().toInt(), Jss = JssB->currentText().toInt();
	int I = Iso->currentIndex();
	QStringList HL;
	aPQR = Js - Jss + 1;
	Tab->setRowCount(Nvs[aPQR][I][Jss]);
	Tab->setColumnCount(Nvss[aPQR][I][Jss]);
	for (vss = 0; vss < Nvss[aPQR][I][Jss]; vss++) HL << "v''=" + QString::number(vss);
	Tab->setHorizontalHeaderLabels(HL);
	HL.clear();
	for (vs = 0; vs < Nvs[aPQR][I][Jss]; vs++) HL << "v'=" + QString::number(vs);
	Tab->setVerticalHeaderLabels(HL);
	for (vs = 0; vs < Nvs[aPQR][I][Jss]; vs++) for (vss = 0; vss < Nvss[aPQR][I][Jss]; vss++)
	{
		if (Tab->item(vs, vss) != 0) 
			Tab->item(vs, vss)->setText(QString::number(Data[aPQR][I][Jss][vs][vss], 'f', 6));
		else Tab->setItem(vs, vss, 
			new QTableWidgetItem(QString::number(Data[aPQR][I][Jss][vs][vss], 'f', 6)));
	}
}

bool FCFTab::readData(QString Filename)
{
	QFile Datei(Filename);
	if (!read(&Datei)) return false;
	QDataStream S(&Datei);
	char *s = new char[7];
	int l, i, I, J, vs, vss;
	bool Err = false;
	if (S.readRawData(s, 6) < 6) return false;
	s[6] = 0;
	if (s != QString("FCFTab")) return false;
	delete[] s;
	S >> l;
	if (l<=0) return false;
	s = new char[l+1];
	if (S.readRawData(s, l) < l) return false;
	s[l] = 0;
	setName(s);
	S >> l;
	if (l<=0) return false;
	s = new char[l+1];
	if (S.readRawData(s, l) < l) return false;
	s[l] = 0;
	Pot1->setText(s);
	S >> l;
	if (l<=0) return false;
	s = new char[l+1];
	if (S.readRawData(s, l) < l) return false;
	s[l] = 0;
	Pot2->setText(s);
	S >> l;
	if (l<=0) return false;
	s = new char[l+1];
	if (S.readRawData(s, l) < l) return false;
	s[l] = 0;
	Date->setText(s);
	DestroyData();
	NIso = new int[3];
	NJ = new int*[3];
	Nvs = new int**[3];
	Nvss = new int**[3];
	Data = new float****[3];
	for (i=0; i<3; i++)
	{
		if (Err) NIso[i] = 0;
		else S >> NIso[i];
		if (NIso[i] < 0) 
		{
			NIso[i] = 0;
			Err = true;
		}
		if (NIso[i] == 0)
		{
			NJ[i] = 0;
			Nvs[i] = Nvss[i] = 0;
			Data[i] = 0;
			continue;
		}
		NJ[i] = new int[NIso[i]];
		Nvs[i] = new int*[NIso[i]];
		Nvss[i] = new int*[NIso[i]];
		Data[i] = new float***[NIso[i]];
		for (I=0; I < NIso[i]; I++)
		{
			if (Err) NJ[i][I] = 0;
			else S >> NJ[i][I];
			if (NJ[i][I] < 0) 
			{
				NJ[i][I] = 0;
				Err = true;
			}
			if (NJ[i][I] == 0)
			{
				Nvs[i] = Nvss[i] = 0;
				Data[i] = 0;
				continue;
			}
			Nvs[i][I] = new int[NJ[i][I]];
			Nvss[i][I] = new int[NJ[i][I]];
			Data[i][I] = new float**[NJ[i][I]];
			for (J=0; J < NJ[i][I]; J++)
			{
				if (Err) Nvs[i][I][J] = Nvss[i][I][J] = 0;
				else S >> Nvs[i][I][J] >> Nvss[i][I][J];
				if (Nvs[i][I][J] < 0 || Nvss[i][I][J] < 0)
				{
					Nvs[i][I][J] = 0;
					Err = true;
				}
				if (Nvs[i][I][J] == 0 || Nvss[i][I][J] == 0)
				{
					Nvss[i][I][J] = Nvs[i][I][J] = 0;
					Data[i][I][J] = 0;
					continue;
				}
				Data[i][I][J] = new float*[Nvs[i][I][J]];
				for (vs=0; vs < Nvs[i][I][J]; vs++)
				{
					Data[i][I][J][vs] = new float[Nvss[i][I][J]];
					for (vss = 0; vss < Nvss[i][I][J]; vss++) S >> Data[i][I][J][vs][vss];
				}
			}
		}
	}
	if (Err) return false;
	return true;
}

void FCFTab::setData(const int PQR, const int Iso, const int J, const int Nvs, const int Nvss, double **** FCF)
{
    if (Nvs == 0 || Nvss == 0)
    {
        setData(PQR, Iso, J, 0, 0, static_cast<double***>(NULL));
        return;
    }
    double*** acFCF = Create(Nvs, Nvss, 2), *ChannelF = (PQR == 1 ? m_ChannelFQ : m_ChannelFPR), dnull = 0.0;
    int Nchannels = (PQR == 1 ? m_NChannelsQ : m_NChannelsPR);
    for (int vs = 0; vs < Nvs; ++vs)
    {
        for (int vss = 0; vss < Nvss; ++vss)
        {
            acFCF[vs][vss][0] = dnull;
            for (int c=0; c < Nchannels; ++c) acFCF[vs][vss][0] += FCF[c][vs][vss][0] * ChannelF[c];
            acFCF[vs][vss][1] = FCF[0][vs][vss][1];
            for (int c=0; c < Nchannels; ++c) delete[] FCF[c][vs][vss];
        }
        for (int c=0; c < Nchannels; ++c) delete[] FCF[c][vs];
    }
    for (int c=0; c < Nchannels; ++c) delete[] FCF[c];
    delete[] FCF;
    setData(PQR, Iso, J, Nvs, Nvss, acFCF);
}

void FCFTab::setData(int PQR, int I, int J, int nvs, int nvss, double*** FCF)
{
	if (nvs == 0 || nvss == 0)
	{
		Nvs[PQR][I][J] = Nvss[PQR][I][J] = 0;
		Data[PQR][I][J] = 0;
		return;
	}
	int vs, vss;
	Data[PQR][I][J] = new float*[nvs];
	for (vs = 0; vs < nvs; vs++)
	{
		Data[PQR][I][J][vs] = new float[nvss];
		for (vss = 0; vss < nvss; vss++) 
		{
			Data[PQR][I][J][vs][vss] = float(FCF[vs][vss][0]);
			delete[] FCF[vs][vss];
		}
		delete[] FCF[vs];
	}
	delete[] FCF;
}

void FCFTab::setData(QString pot1, QString pot2, int* nIso, int** nJ, int*** nvs, int*** nvss,
					 float***** data)
{
	DestroyData();
	Pot1->setText(pot1);
	Pot2->setText(pot2);
	Date->setText(QDateTime::currentDateTime().toString("dd.MM.yyyy hh:mm"));
	NIso = nIso;
	NJ = nJ;
	Nvs = nvs;
	Nvss = nvss;
	Data = data;
	Changed();
}

void FCFTab::setStates(ElState* State1, ElState* State2)
{
	St1 = State1;
	St2 = State2;
}

void FCFTab::setTransition(Transition* T)
{
	Trans = T;
	St2 = T->getLowerState();
	St1 = T->getUpperState();
}

void FCFTab::shiftPar(double *&E1, double *&E2, double ***&WF1, double ***&WF2, int *&ri1,
					  int *&ri2, int *&ra1, int *&ra2, int *&P1, int *&P2, double *&F1, 
					  double *&F2, double *&N1, double *&N2, int &Nv1, int &Nv2)
{
	E1 = E2;
	WF1 = WF2;
	ri1 = ri2;
	ra1 = ra2;
	P1 = P2;
	F1 = F2;
	N1 = N2;
	Nv1 = Nv2;
	E2 = 0;
	WF2 = 0;
	ri2 = 0;
	ra2 = 0;
	P2 = 0;
	F2 = 0;
	N2 = 0;
	Nv2 = 0;
}

bool FCFTab::writeData(QString Filename)
{
	QFile Datei(Filename);
	if (!write(&Datei)) return false;
	QDataStream S(&Datei);
	QString T = getName();
	int l = T.length(), i, I, J, vs, vss;
	if (S.writeRawData("FCFTab", 6) < 6) return false;
	S << l;
	if (S.writeRawData(T.toLatin1().constData(), l) < l) return false;
	S << (l = (T = Pot1->text()).length());
	if (S.writeRawData(T.toLatin1().constData(), l) < l) return false;
	S << (l = (T = Pot2->text()).length());
	if (S.writeRawData(T.toLatin1().constData(), l) < l) return false;
	S << (l = (T = Date->text()).length());
	if (S.writeRawData(T.toLatin1().constData(), l) < l) return false;
	for (i=0; i<3; i++)
	{
		S << NIso[i];
		for (I=0; I < NIso[i]; I++)
		{
			S << NJ[i][I];
			for (J=0; J < NJ[i][I]; J++)
			{
				S << Nvs[i][I][J] << Nvss[i][I][J];
				for (vs = 0; vs < Nvs[i][I][J]; vs++) for (vss = 0; vss < Nvss[i][I][J]; vss++)
					S << Data[i][I][J][vs][vss];
			}
		}
	}
	if (Datei.error() != QFile::NoError) return false;
	return true;
}


