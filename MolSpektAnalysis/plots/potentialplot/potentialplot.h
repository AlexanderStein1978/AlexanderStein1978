//
// C++ Interface: potentialplot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef POTENTIALPLOT_H
#define POTENTIALPLOT_H

#include "DiagWindow.h"

class Potential;
class MainWindow;
struct SplinePoint;


class PotentialPlot : public DiagWindow
{
	Q_OBJECT

public:
    PotentialPlot(Potential *potential = 0, MainWindow *MW = 0);
    ~PotentialPlot();
	void plotPotential(Potential *potential = 0);
	void SetBoundaries();
    void PSpektrum(QPainter &P, const QRect & A, bool PrintFN );
	void ResetZoom();
	void clearHistory();
	void addPotential(Potential *pot);
	void setShowHistory(bool H);
	bool getShowHistory();
	void PotSnapShot();
	bool getShowDiagFuncs();
	void setShowDiagFuncs(bool show);
	
public slots:
	void PotentialChanged();
	
private slots:
	void addPoint();
	void removePoint();
	
private:
    void SetPoints();
    void MovePoint();
    void MovePoint(int i_pointIndex, double i_newX, double i_newY);
    void ShowPopupMenu(const QPoint &i_point);
    void HandleHistoryWhileMoving();
    void HandleHistoryAfterMoving();

	Potential *pot, **potCopies;
    int NCopies, copyP;
    bool History, showDiagFuncs;
    int nDFP;
	QMenu *pointMenu, *pictureMenu;
	QAction *removePointAct, *addPointAct;
    double *WFS, *WWFS, DFRi, DFh;
};

#endif
