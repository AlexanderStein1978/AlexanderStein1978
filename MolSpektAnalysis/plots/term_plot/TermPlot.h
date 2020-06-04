//
// C++ Interface: TermPlot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TERMPLOT_H
#define TERMPLOT_H

#include <qwidget.h>
#include <qlineedit.h>
#include "tools.h"
#include "DiagWindow.h"

class Tabelle;
class Molecule;
class ElState;
class LineTable;


class TermPlot : public DiagWindow
{
    Q_OBJECT
public:    
    TermPlot(MainWindow *MW = 0, Molecule *Mol = 0, ElState *State = 0);
    ~TermPlot();
	void SetTESources();
public slots:
	void Zoom();
	//void Mark();
	//void SMark(int x1, int y1, int x2, int y2);
    void Updatevs();
    void UpdateLines();
	void UpdateMarked();
    //void PlotDunham(int KMax, int *LMax, double **Par, double **Korr);
    void SetBoundaries();
    void PSpektrum(QPainter &P, const QRect & A, bool PrintFN );
	void PictureClicked(QPoint *P);
signals:
    void vsUpdated();
private slots:
	void IsoChanged(int nIso);	
private:
    void drawP(int x, int y, QColor C, QPainter &p);
    int slx1, sly1, slx2, sly2, *vs, *vss, *Jss, *MIso, *LIso;
	int iso, nlC, nlI, nlv, nlJ, *LJtJp1, *MJtJp1, AnzahlLinien, AnzahlMarked;
	double ****ELU, *LTE, *MTE;
	bool SMarked, *LMarked;
	Molecule *mol;
	ElState *state, *lstate;
	LineTable **LLT;
};
	    
#endif
    
