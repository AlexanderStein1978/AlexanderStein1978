//
// C++ Interface: AddSpectrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef ADDSPEKTRUM_H
#define ADDSPEKTRUM_H


#include "DiagWindow.h"

#include <QToolButton>
#include <QMenu>

class AddDialog;
class Spektrum;
class ElState;
class TermTable;
class FitData;
class DoubleAddSpectrum;
class ResidualFit;

struct MarkedPeak;
struct TermEnergy;


class AddSpectrum : public DiagWindow
{
	Q_OBJECT
	
public:
	AddSpectrum(MainWindow *MW, AddDialog *D);
	~AddSpectrum();
	void setData(double **nData, int NJ, int JStart, int JStep, int NE, Spektrum *Source, 
				 ElState *State, int Comp, int Iso, int v, TermTable *TT, double Ri, 
                 double Ra, double Res, double *UEnergy, int *vAssign, int *cAssign, double *i_IntensityF);
	void getMarked(int &N, MarkedPeak *&MarkedLevels);
	void maximizeDAdd(double Diff);
	int AutoUpdateLevels(bool fullyAutomatic);
	
	inline int getSpectrumType()
	{
		return type;
	}
	
public slots:
	virtual void PSpektrum(QPainter &P, const QRect &A, bool PrintFN);
	virtual void KeyPressEvent(QKeyEvent *K);
	void acceptAutoUpdateResult();
	void rejectAutoUpdateResult();
	
private slots:
	void mark(QRect *R);
	void removeLevels();
	void updateLevels();
	void hideLevels();
	void showLevels();
	void markLevels();
	void setEnergyReference();
	void ImageClicked(double dJ, double E);
	void AssignmentsAccepted(FitData *CFD);
	void stateChanged();
	void drawDAddSpectrum(double dRes, bool newState, int SplineNum = -1, int JMin = -1, int JMax = -1);
	void fillTEStruct(TermEnergy &TE, int &FDatNr, int JPos, int EPos);
	void addLevel(TermEnergy &TE, int &FDatNr, int JPos, int EPos);
    void contrastChanged(int i_c);
	
	inline void autoUpdateLevels()
	{
		AutoUpdateLevels(false);
	}
	
	inline void showPMenu(QPoint* P)
	{
		if (!ZoomB->isChecked())
		{
			if (NMBAUL == -1) PMenu->popup(*P);
			else AutoUpdateMenu->popup(*P);
		}
	}
	
private:
	void calcEJ(double &rE, int &rJ, int J, int P);
	int getPeakTop(int JPos, int StartPos);
	
	FitData *FitD, **PFitD;
	MarkedPeak *marked;
	int NMarked, type, Comp, TDC, v, Iso, NP, NJ, JStart, JStep, NC, NMBAUL;
    int NPert, *NCP, *NvP, *NJP, *EShift, NOC, NOIso, NOJ, NOv, *vAss, *CAss, m_contrast;
    double **Data, Ri, Ra, Res, *****PTData, ****OTData, *UEn, *m_IntensF;
	bool **UCP, m_ownResidualFit;
	TermTable *TT, *OffTT, *PTT[50];
	ElState *State;
	Spektrum *Source;
    AddDialog *Dialog;
	DoubleAddSpectrum *DAddSpect;
	ResidualFit *ResFit;
	QMenu *PMenu, *AutoUpdateMenu;
	QAction *removeLevelsAct, *updateLevelsAct, *hideLevelsAct, *showLevelsAct;
	QAction *markLevelsAct, *autoUpdateAct;
};

#endif
