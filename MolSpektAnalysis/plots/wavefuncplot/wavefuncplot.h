//
// C++ Interface: wavefuncplot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef WAVEFUNCPLOT_H
#define WAVEFUNCPLOT_H

#include "DiagWindow.h"


class CoupledSineWaveFunc;
class MainWindow;
class Potential;
class Molecule;

class QPaintDevice;
class QRect;
class QColor;

class WaveFuncPlot : public DiagWindow
{
	Q_OBJECT
			
	public:
        WaveFuncPlot(int NumWFPoints, MainWindow *MW = 0, Potential *Pot = 0, int Iso = -1, int J = 0, int v = -1);
	    ~WaveFuncPlot();
		
	private slots:
		void PotChanged(QString Name);
		void IsoChanged(int I);
		void JChanged();
		void vChanged(int v);
		
	private:
		void PSpektrum(QPainter &P, const QRect &A, bool PrintFN);
		void Free();
		
        int *P, *ri, *ra, Nv, NChannels, NumWFPoints;
		double WaveFuncSF, ***WF, *E, *F, *N, h, Ra, Ro, PMin, PMax;
		QColor *PotC, *RotPotC, *LevelC, *WaveC;
		Molecule *Mol;
		Potential **Pot;
		CoupledSineWaveFunc *SWF;
};

#endif
